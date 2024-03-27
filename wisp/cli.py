import argparse
import os

from . import __version__
from .contexts import ContextManager


def setup_cli_interface(context_manager: ContextManager) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Update parameters for WISP.")
    attributes = context_manager.get()  # Assuming this returns a dict of attributes

    parser.add_argument("pdb_path", type=str, help="Path to PDB file to analyze")
    for attr, value in attributes.items():
        # Use the type of the current value to infer the expected type
        # This is a simplified approach; you might need custom handling for complex types
        arg_type = type(value) if value is not None else str
        if attr in ["source_residues", "sink_residues"]:
            parser.add_argument(
                f"--{attr}",
                nargs="+",
                required=True,
            )
        else:
            parser.add_argument(
                f"--{attr}", type=arg_type, help=f"Set {attr} (current value: {value})"
            )

    return parser


def setup_output_dir(output_dir: str) -> None:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


def run_cli() -> ContextManager:
    context_manager = ContextManager()
    parser = setup_cli_interface(context_manager)
    args = parser.parse_args()

    # Convert argparse Namespace to dictionary, excluding None values
    updates = {k: v for k, v in vars(args).items() if v is not None}
    context_manager.update(updates)

    setup_output_dir(context_manager.output_dir)
    return context_manager


if __name__ == "__main__":
    run_cli()
