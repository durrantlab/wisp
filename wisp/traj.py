import multiprocessing as mp

from loguru import logger

from .structure import Molecule


# --- first, to generate a covariance matrix ---#
class multi_threading_to_collect_data_from_frames:
    """Launch PDB-frame processing on multiple processors"""

    combined_results = None

    def __init__(self, inputs, num_processors: int | None = None):
        """
        Args:
            inputs: the data to be processed, in a list
            num_processors: the number of processors to use to process this data,
                an integer

        """
        self.results = []

        # First, we determine the number of available cores.
        if num_processors is None:
            num_processors = mp.cpu_count()
            logger.debug("Setting the number of cores to ", num_processors)

        # reduce the number of processors if too many have been specified
        if len(inputs) < num_processors:
            logger.debug("Number of cores is higher than number of inputs.")
            num_processors = len(inputs)
            if num_processors == 0:
                num_processors = 1
            logger.debug("Setting number of cores to ", num_processors)

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
            threads.append(collect_data_from_frames())
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
        total_summed_coordinates = None
        dictionary_of_node_lists = {}
        for _ in threads:
            chunk = results_queue.get()

            if total_summed_coordinates is None:
                total_summed_coordinates = chunk[0]
            else:
                total_summed_coordinates = total_summed_coordinates + chunk[0]

            for key in chunk[1].keys():
                try:
                    dictionary_of_node_lists[key].extend(chunk[1][key])
                except Exception:
                    dictionary_of_node_lists[key] = chunk[1][key]

        self.combined_results = (total_summed_coordinates, dictionary_of_node_lists)


class collect_data_from_frames:
    """PDB-frame data processing on a single processor"""

    summed_coordinates = None
    nodes = {}

    def runit(self, running, mutex, results_queue, items):
        """
        Args:
            running: a mp.Value() object
            mutex: a mp.Lock() object
            results_queue: where the results will be stored [mp.Queue()]
            items: the data to be processed, in a list
        """
        for item in items:
            self.value_func(item)  # , results_queue)
        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put((self.summed_coordinates, self.nodes))

    def value_func(
        self, params_and_res_keys_and_pdb_lines_and_res_maps
    ):  # , results_queue): # so overwriting this function
        """Process a single PDB frame: identify the relevant nodes

        Args:
            params_and_res_keys_and_pdb_lines_and_res_maps: a tuple containing required
                information.
                The first item contains user-defined parameters (a UserInput object)
                The second item is a list containing string representations of each residue ("CHAIN_RESNAME_RESID")
                The third item is a list of strings representing the PDB frame to be processed, where each string
                    contains a PDB ATOM or HETATM entry
                The fourth item is a dictionary that maps residue string identifiers ("CHAIN_RESNAME_RESID") to a list
                    of the indices of the atoms that correspond to that residue
        """

        params = params_and_res_keys_and_pdb_lines_and_res_maps[
            0
        ]  # user-defined parameters
        pdb_lines = params_and_res_keys_and_pdb_lines_and_res_maps[
            1
        ]  # make sure this is not empty

        # now load the frame into its own Molecule object
        pdb = Molecule()
        pdb.load_pdb_from_list(pdb_lines)

        if self.summed_coordinates is None:
            self.summed_coordinates = pdb.coordinates
        else:
            self.summed_coordinates = self.summed_coordinates + pdb.coordinates

        pdb.map_atoms_to_residues()
        pdb.map_nodes_to_residues(params["node_definition"])

        for index, residue_iden in enumerate(pdb.residue_identifiers_in_order):
            try:
                self.nodes[residue_iden].append(pdb.nodes[index])
            except Exception:
                self.nodes[residue_iden] = [pdb.nodes[index]]
