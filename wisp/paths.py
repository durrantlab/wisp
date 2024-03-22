import copy
import multiprocessing as mp
import sys
import time
from collections.abc import Collection

import networkx as nx
import numpy as np
from loguru import logger


def get_log_n_paths(graph, cutoff_length):
    # Calculate the average branching factor
    total_edges = graph.number_of_edges()
    total_nodes = graph.number_of_nodes()
    avg_branching_factor = total_edges / total_nodes if total_nodes else 0

    # Use logarithms to avoid overflow
    # Check if avg_branching_factor is greater than 1 to avoid log(0) or negative values
    if avg_branching_factor > 1:
        log_estimated_paths = cutoff_length * np.log(avg_branching_factor)
    else:
        # If the avg_branching_factor is 1 or less, the growth is linear or
        # non-existent, not exponential
        log_estimated_paths = 0

    logger.debug(f"log(estimated_n_paths) = {log_estimated_paths}")

    return log_estimated_paths


class multi_threading_find_paths:
    """Launches path finding on multiple processors"""

    results = []

    def __init__(self, inputs: Collection, num_processors: int | None = None):
        """
        Args:
            inputs: the data to be processed, in a list
            num_processors: the number of processors to use to process this data.
        """

        self.results = []

        # First, we determine the number of available cores.
        if num_processors is None:
            num_processors = mp.cpu_count()
        # reduce the number of processors if too many have been specified
        if len(inputs) < num_processors:
            logger.debug("Number of cores is higher than number of inputs.")
            num_processors = len(inputs)
            if num_processors == 0:
                num_processors = 1
        logger.debug(f"Setting the number of cores to {num_processors}")

        # now, divide the inputs into the appropriate number of processors
        inputs_divided = {t: [] for t in range(num_processors)}
        for t in range(0, len(inputs), num_processors):
            for t2 in range(num_processors):
                index = t + t2
                if index < len(inputs):
                    inputs_divided[t2].append(inputs[index])

        # now, run each division on its own processor
        running = mp.Value("i", num_processors)
        mutex = mp.Lock()

        arrays = []
        threads = []
        for _ in range(num_processors):
            threads.append(find_paths())
            arrays.append(mp.Array("i", [0, 1]))

        results_queue = mp.Queue()  # to keep track of the results

        processes = []
        for i in range(num_processors):
            p = mp.Process(
                target=threads[i].runit,
                args=(running, mutex, results_queue, inputs_divided[i]),
            )
            p.start()
            processes.append(p)

        while running.value > 0:
            continue  # wait for everything to finish

        # compile all results
        for _ in threads:
            chunk = results_queue.get()
            for chun in chunk:
                self.results.extend(chun)


class find_paths:  # other, more specific classes with inherit this one
    """Path-finding data processing on a single processor"""

    results = []

    def runit(self, running, mutex, results_queue, items):
        """Path-finding data processing on a single processor.

        Args:
            running: a mp.Value() object
            mutex: a mp.Lock() object
            results_queue: where the results will be stored [mp.Queue()]
            items: the data to be processed, in a list
        """

        for item in items:
            self.value_func(item, results_queue)
        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put(self.results)

    def value_func(
        self, item, results_queue
    ):  # this is the function that changes through inheritance
        """Process a single path-finding "branch"

        Args:
            item: a tuple containing required information.
                The first is a numpy array containing a single float, the path-length cutoff
                The second is an index corresponding to the ultimate path sink
                The third is a nx.Graph object describing the connectivity of the different nodes
                The fourth is a list corresponding to a path. The first item is the length of the path (float).
                    The remaining items are the indices of the nodes in the path (int).
            results_queue: where the results will be stored [mp.Queue()]
        """

        cutoff = item[0]
        sink = item[1]
        G = item[2]

        paths_growing_out_from_source = [item[3]]
        full_paths_from_start_to_sink = []

        while paths_growing_out_from_source:
            self.expand_growing_paths_one_step(
                paths_growing_out_from_source,
                full_paths_from_start_to_sink,
                cutoff,
                sink,
                G,
            )

        # here save the results for later compilation
        self.results.append(full_paths_from_start_to_sink)

    def expand_growing_paths_one_step(
        self,
        paths_growing_out_from_source,
        full_paths_from_start_to_sink,
        cutoff,
        sink,
        G,
    ):
        """Expand the paths growing out from the source to the sink by one step
           (to the neighbors of the terminal node) of the expanding paths

        Args:
            paths_growing_out_from_source: a list of paths, where each path is
                represented by a list. The first item in each path is the length of
                the path (float). The remaining items are the indices of the nodes
                in the path (int).
            full_paths_from_start_to_sink: a growing list of identified paths that
                connect the source and the sink, where each path is formatted as above.
            cutoff: a numpy array containing a single element (float), the length
                cutoff. Paths with lengths greater than the cutoff will be ignored.
            sink: the index of the sink (int)
            G: a nx.Graph object describing the connectivity of the different nodes
        """

        for i, path_growing_out_from_source in enumerate(paths_growing_out_from_source):
            if path_growing_out_from_source[0] > cutoff:
                # Because if the path is already greater than the cutoff, no
                # use continuing to branch out, since subsequent branhes will
                # be longer.
                paths_growing_out_from_source.pop(i)
                break
            elif path_growing_out_from_source[-1] == sink:
                # so the sink has been reached
                full_paths_from_start_to_sink.append(path_growing_out_from_source)
                paths_growing_out_from_source.pop(i)
                break
            elif path_growing_out_from_source[-1] != sink:
                # sink not yet reached, but paths still short enough. So add
                # new paths, same as old, but with neighboring element
                # appended.
                node_neighbors = list(G.neighbors(path_growing_out_from_source[-1]))
                for j, node_neighbor in enumerate(node_neighbors):
                    if not node_neighbor in path_growing_out_from_source:
                        temp = path_growing_out_from_source[:]
                        temp.append(node_neighbor)
                        temp[0] = temp[0] + G.edges[temp[-2], temp[-1]]["weight"]
                        paths_growing_out_from_source.insert((i + j + 1), temp)
                paths_growing_out_from_source.pop(i)
                break
            else:
                logger.critical("SOMETHING IS WRONG")


class GetPaths:
    """Get the paths from a list of sources to a list of sinks"""

    def __init__(
        self, corr_matrix, srcs, snks, params, residue_keys, n_paths_max=1000000
    ):
        """Identify paths that link the source and the sink and order them by their
        lengths.

        Args:
            corr_matrix: a np.array, the calculated correlation matrix
            srcs: a list of ints, the indices of the sources for path finding
            snks: a list of ints, the indices of the sinks for path finding
            params: the user-specified command-line parameters, a UserInput object
            residue_keys: a list containing string representations of each residue
            n_paths_cutoff: Specifies the maximum number of paths to proceed. If we
                estimate the number of paths to be larger than this number, we will
                terminate the calculation.
        """
        if "n_paths_max" in params.keys():
            n_paths_max = params["n_paths_max"]

        # populate graph nodes and weighted edges
        G = nx.Graph(incoming_graph_data=corr_matrix)

        # first calculate length of shortest path between any source and sink
        logger.info("Calculating paths...", params["logfile"])
        logger.info(
            "Calculating the shortest path between any of the specified sources and any of the specified sinks...",
            params["logfile"],
        )
        shortest_length, shortest_path = self.get_shortest_path_length(
            corr_matrix, srcs, snks, G
        )
        logger.info(
            f"The shortest path has length {str(shortest_length)}",
            params["logfile"],
        )

        path = [shortest_length]
        path.extend(shortest_path)
        pths = [
            path
        ]  # need to create this initial path in case only one path is requrested

        cutoff = shortest_length

        # Check for comb explosion
        log_n_paths = get_log_n_paths(G, cutoff)
        if log_n_paths > np.log(n_paths_max):
            logger.error(f"Estimated number of paths is greater than {n_paths_max}")
            logger.error("Please increase n_paths_max to proceed.")
            logger.error("Terminating calculation.")
            sys.exit(1)

        cutoff_yields_max_num_paths_below_target = 0
        cutoff_yields_min_num_paths_above_target = 1000000.0

        # first step, keep incrementing a little until you have more than the desired number of paths
        logger.info(
            "Identifying the cutoff required to produce "
            + str(params["desired_number_of_paths"])
            + " paths...",
            params["logfile"],
        )
        num_paths = 1
        while num_paths < params["desired_number_of_paths"]:
            logger.info(f"Testing the cutoff {str(cutoff)}...", params["logfile"])
            cutoff_in_array = np.array([cutoff], np.float64)
            pths = self.remove_redundant_paths(
                self.get_paths_between_multiple_endpoints(
                    cutoff_in_array, corr_matrix, srcs, snks, G, params
                )
            )
            num_paths = len(pths)

            logger.info(
                f"The cutoff {str(cutoff)} produces {num_paths} paths...",
                params["logfile"],
            )

            if (
                num_paths < params["desired_number_of_paths"]
                and cutoff > cutoff_yields_max_num_paths_below_target
            ):
                cutoff_yields_max_num_paths_below_target = cutoff
            if (
                num_paths > params["desired_number_of_paths"]
                and cutoff < cutoff_yields_min_num_paths_above_target
            ):
                cutoff_yields_min_num_paths_above_target = cutoff

            # Original code adds .1 each time... but this may be to fast for
            # some systems and very slow for others... lets try increasing by
            # a percentage of the minimum path length instead... ideally this
            # could be an input parameter in the future.
            cutoff = cutoff + shortest_length * 0.1

        pths = self.remove_redundant_paths(pths)

        pths.sort()  # sort the paths by length

        if (
            num_paths != params["desired_number_of_paths"]
        ):  # so further refinement is needed
            pths = pths[: params["desired_number_of_paths"]]
            logger.info(
                "Keeping the first "
                + str(params["desired_number_of_paths"])
                + " of these paths...",
                params["logfile"],
            )

        self.paths_description = ""

        self.paths_description = (
            self.paths_description + "\n# Output identified paths" + "\n"
        )
        index = 1

        if params["simply_formatted_paths_filename"] != "":
            simp = open(params["simply_formatted_paths_filename"], "w")
        for path in pths:
            self.paths_description = (
                f"{self.paths_description}Path {str(index)}:" + "\n"
            )
            self.paths_description = (
                f"{self.paths_description}   Length: {str(path[0])}" + "\n"
            )
            self.paths_description = (
                f"{self.paths_description}   Nodes: "
                + " - ".join([residue_keys[item] for item in path[1:]])
                + "\n"
            )
            if params["simply_formatted_paths_filename"] != "":
                simp.write(" ".join([str(item) for item in path]) + "\n")
            index = index + 1
        if params["simply_formatted_paths_filename"] != "":
            simp.close()

        self.paths = pths

    def remove_redundant_paths(self, pths):
        """Removes redundant paths

        Args:
            pths: a list of paths

        Returns:
            A list of paths with the redundant ones eliminated.
        """

        if len(pths) == 1:
            # no reason to check if there's only one
            return pths

        for indx1 in range(len(pths) - 1):
            path1 = pths[indx1]
            if path1 is not None:
                for indx2 in range(indx1 + 1, len(pths)):
                    path2 = pths[indx2]
                    if path2 is not None and len(path1) == len(
                        path2
                    ):  # paths are the same length
                        pth1 = copy.deepcopy(path1[1:])
                        pth2 = copy.deepcopy(path2[1:])

                        if pth1[0] < pth1[-1]:
                            pth1.reverse()
                        if pth2[0] < pth2[-1]:
                            pth2.reverse()

                        if pth1 == pth2:
                            pths[indx2] = None

        while None in pths:
            pths.remove(None)

        return pths

    def get_shortest_path_length(
        self, corr_matrix, srcs, snks, G
    ):  # where sources and sinks are lists
        """Identify the length of the shortest path connecting any of the sources and any of the sinks

        Args:
            corr_matrix: a np.array, the calculated correlation matrix
            srcs: a list of ints, the indices of the sources for path finding
            snks: a list of ints, the indices of the sinks for path finding
            G: a nx.Graph object describing the connectivity of the different nodes

        Returns:
            a float, the length of the shortest path, and a list of ints corresponding
            to the nodes of the shortest path.
        """

        shortest_length = 99999999.999
        shortest_path = []

        for source in srcs:
            for sink in snks:
                if source != sink:  # important to avoid this situation
                    short_path = nx.dijkstra_path(G, source, sink, weight="weight")
                    length = self.get_length_of_path(short_path, corr_matrix)
                    if length < shortest_length:
                        shortest_length = length
                        shortest_path = short_path
        return shortest_length, shortest_path

    def get_length_of_path(self, path, corr_matrix):
        """Calculate the length of a path

        Args:
            path: a list of ints, the indices of the path
            corr_matrix: a np.array, the calculated correlation matrix

        Returns:
            a float, the length of the path
        """

        length = 0.0
        for t in range(len(path) - 1):
            length = length + corr_matrix[path[t], path[t + 1]]
        return length

    def get_paths_between_multiple_endpoints(
        self, cutoff, corr_matrix, srcs, snks, G, params
    ):  # where sources and sinks are lists
        """Get paths between sinks and sources

        Args:
            cutoff: a np.array containing a single float, the cutoffspecifying the maximum permissible path length
            corr_matrix: a np.array, the calculated correlation matrix
            srcs: a list of ints, the indices of the sources for path finding
            snks: a list of ints, the indices of the sinks for path finding
            G: a nx.Graph object describing the connectivity of the different nodes
            params: the user-specified command-line parameters, a UserInput object

        Returns:
            a list of paths, where each path is represented by a list. The first item in each path is the length
            of the path (float). The remaining items are the indices of the nodes in the path (int).
        """

        pths = []
        for source in srcs:
            for sink in snks:
                if source != sink:  # avoid this situation
                    pths.extend(
                        self.get_paths_fixed_endpoints(
                            cutoff, corr_matrix, source, sink, G, params
                        )
                    )
        return pths

    def get_paths_fixed_endpoints(self, cutoff, corr_matrix, source, sink, G, params):
        """Get paths between a single sink and a single source

        Args:
            cutoff: a np.array containing a single float, the cutoff specifying the
                maximum permissible path length
            corr_matrix: a np.array, the calculated correlation matrix
            source: the index of the source for path finding
            sink: the index of the sink for path finding
            G: a nx.Graph object describing the connectivity of the different nodes
            params: the user-specified command-line parameters, a UserInput object

        Returns:
            a list of paths, where each path is represented by a list. The first item
            in each path is the length of the path (float). The remaining items are
            the indices of the nodes in the path (int).
        """

        if source == sink:
            return []

        source_lengths, source_paths = nx.single_source_dijkstra(
            G, source, target=None, cutoff=None, weight="weight"
        )
        sink_lengths, sink_paths = nx.single_source_dijkstra(
            G, sink, target=None, cutoff=None, weight="weight"
        )

        so_l = [source_lengths[key] for key in source_lengths.keys()]
        so_p = [source_paths[key] for key in source_paths.keys()]
        si_l = [sink_lengths[key] for key in sink_lengths.keys()]
        si_p = [sink_paths[key] for key in sink_paths.keys()]

        check_list_1 = []
        check_list_2 = []
        for i in range(len(so_l)):
            check_list_1.extend([so_p[i][-1]])
            check_list_2.extend([si_p[i][-1]])

        node_list = []
        dijkstra_list = []
        upper_minimum_length = 0
        if not set(check_list_1).difference(check_list_2):
            for i, _ in enumerate(so_l):
                if so_l[i] + si_l[i] <= cutoff:
                    node_list.extend(so_p[i][:])
                    node_list.extend(si_p[i][:])
                    si_pReversed = si_p[i][:]
                    si_pReversed.reverse()
                    temp_path = so_p[i][:] + si_pReversed[1:]
                    temp_length = so_l[i] + si_l[i]
                    dijkstra_list.append(temp_path)
                    if (so_l[i] + si_l[i]) > upper_minimum_length:
                        upper_minimum_length = temp_length
        else:
            logger.critical("paths do not match up")

        unique_nodes = list(set(node_list))
        unique_nodes.sort()

        node_length = len(unique_nodes)
        new_matrix = np.zeros((len(corr_matrix), len(corr_matrix)))

        for i in range(node_length):
            for j in range(node_length):
                new_matrix[unique_nodes[i]][unique_nodes[j]] = corr_matrix[
                    unique_nodes[i]
                ][unique_nodes[j]]

        corr_matrix = new_matrix
        G = nx.Graph(incoming_graph_data=corr_matrix, labels=unique_nodes)

        length = 0.0
        paths_growing_out_from_source = [[length, source]]
        full_paths_from_start_to_sink = []

        # This is essentially this list-addition replacement for a recursive
        # algorithm you've envisioned.
        # To parallelize, just get the first N branches, and send them off to each node.
        # Rest of branches filled out in separate processes.

        find_paths_object = find_paths()
        if params["number_processors"] == 1:
            while paths_growing_out_from_source:
                find_paths_object.expand_growing_paths_one_step(
                    paths_growing_out_from_source,
                    full_paths_from_start_to_sink,
                    cutoff,
                    sink,
                    G,
                )
        else:
            # just get some of the initial paths on a single processor
            logger.info(
                "Starting serial portion of path-finding algorithm (will run for "
                + str(params["seconds_to_wait_before_parallelizing_path_finding"])
                + " seconds)...",
                params["logfile"],
            )
            atime = time.time()
            while (
                paths_growing_out_from_source
                and time.time() - atime
                < params["seconds_to_wait_before_parallelizing_path_finding"]
            ):
                find_paths_object.expand_growing_paths_one_step(
                    paths_growing_out_from_source,
                    full_paths_from_start_to_sink,
                    cutoff,
                    sink,
                    G,
                )

            # ok, so having generated just a first few, divy up those among multiple processors
            if paths_growing_out_from_source:  # in case you've already finished
                logger.info(
                    "Starting parallel portion of path-finding algorithm running on "
                    + str(params["number_processors"])
                    + " processors...",
                    params["logfile"],
                )
                paths_growing_out_from_source = [
                    (cutoff, sink, G, path) for path in paths_growing_out_from_source
                ]
                additional_full_paths_from_start_to_sink = multi_threading_find_paths(
                    paths_growing_out_from_source, params["number_processors"]
                )
                full_paths_from_start_to_sink.extend(
                    additional_full_paths_from_start_to_sink.results
                )
            else:
                logger.info(
                    "(All paths found during serial path finding; parallelization not required)",
                    params["logfile"],
                )

        full_paths_from_start_to_sink.sort()

        pths = []

        for full_path_from_start_to_sink in full_paths_from_start_to_sink:
            pths.append(full_path_from_start_to_sink)

        return pths
