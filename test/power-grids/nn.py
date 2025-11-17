#!/usr/bin/env python3
"""Power grid circuit generator."""

import random
import argparse
from tqdm import tqdm

# Voltage values
VDD = 1.8
GND = 0.0

# Component values
CAPACITOR_VALUE = 5.37096111111111e-11  # Capacitor value in Farads
RESISTOR_VALUE = 2.50e-01  # Grid resistor value in Ohms
SOURCE_RESISTOR_VALUE = 2.5e-02  # Source resistor value in Ohms
VIA_RESISTANCE = 0.0  # Via resistance in Ohms (currently 0)

# Current source parameters
CURRENT_SOURCE_VALUE = 2.5099e-05  # Current source value in Amperes
PULSE_LOW = CURRENT_SOURCE_VALUE  # Pulse low value
PULSE_HIGH = 0.0627475  # Pulse high value
PULSE_DELAY = 1.0e-9  # Pulse delay
PULSE_RISE_TIME = 1e-19  # Pulse rise time
PULSE_FALL_TIME = 1e-10  # Pulse fall time
PULSE_ON_TIME = 1e-19  # Pulse on time
PULSE_PERIOD = 3e-09  # Pulse period

# Simulation parameters
TRAN_STEP = 1e-12  # Transient analysis step time
TRAN_STOP = 1e-8  # Transient analysis stop time


def write_via(out_print, l1, l2, x, y, via_resistance, digits):
    """Write a via (vertical connection) between two layers."""
    out_print(
        f"VL{l1}P{x:0{digits}}{y:0{digits}}L{l2}P{x:0{digits}}{y:0{digits}} nL{l1}P{x:0{digits}}{y:0{digits}} nL{l2}P{x:0{digits}}{y:0{digits}} {via_resistance}"
    )


def calculate_default_vias(nx, ny):
    """Calculate default number of vias as 20% of grid size."""
    return int(0.2 * nx * ny)


def get_default_outfile(args):
    """Generate default output filename based on grid parameters."""
    return f"grid_{args.nx}x{args.ny}_l{args.nl}_v{args.num_vias}.sp"


def main():
    parser = argparse.ArgumentParser(
        description="Generate a power grid netlist.", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("outfile", nargs="?", help="Output file path (default: grid_NXxNY_lNL_vNVIAS.sp)")
    parser.add_argument("--nx", type=int, default=4, help="Grid size in X dimension")
    parser.add_argument("--ny", type=int, help="Grid size in Y dimension (default: nx)")
    parser.add_argument("--nl", type=int, default=2, help="Number of layers")
    parser.add_argument("--num-vias", type=int, help="Number of vias to add (default: 20%% of grid size)")

    # Component value arguments
    parser.add_argument("--capacitor-value", type=float, default=CAPACITOR_VALUE, help=f"Capacitor value in Farads")
    parser.add_argument("--resistor-value", type=float, default=RESISTOR_VALUE, help=f"Grid resistor value in Ohms")
    parser.add_argument(
        "--source-resistor-value", type=float, default=SOURCE_RESISTOR_VALUE, help=f"Source resistor value in Ohms"
    )
    parser.add_argument("--via-resistance", type=float, default=VIA_RESISTANCE, help=f"Via resistance in Ohms")
    parser.add_argument(
        "--current-source-value",
        type=float,
        default=CURRENT_SOURCE_VALUE,
        help=f"Current source value in Amperes",
    )
    parser.add_argument(
        "--pulse-high-value",
        type=float,
        default=PULSE_HIGH,
        help=f"Pulse high amplitude value",
    )

    # Pulse timing parameters (grouped together)
    parser.add_argument(
        "--pulse-timing",
        type=float,
        nargs=5,
        default=[PULSE_DELAY, PULSE_RISE_TIME, PULSE_FALL_TIME, PULSE_ON_TIME, PULSE_PERIOD],
        metavar=("DELAY", "RISE_TIME", "FALL_TIME", "ON_TIME", "PERIOD"),
        help=f"Pulse timing parameters: delay rise_time fall_time on_time period ",
    )
    # Simulation parameters (grouped together)
    parser.add_argument(
        "--sim-params",
        type=float,
        nargs=2,
        default=[TRAN_STEP, TRAN_STOP],
        metavar=("STEP", "STOP"),
        help=f"Simulation parameters: step_time stop_time",
    )

    args = parser.parse_args()

    # Assign component values from arguments
    capacitor_value = args.capacitor_value
    resistor_value = args.resistor_value
    source_resistor_value = args.source_resistor_value
    via_resistance = args.via_resistance
    current_source_value = args.current_source_value
    pulse_high_value = args.pulse_high_value

    # Assign pulse timing parameters from arguments
    pulse_low = current_source_value  # Low value equals current source value
    pulse_delay, pulse_rise_time, pulse_fall_time, pulse_on_time, pulse_period = args.pulse_timing

    # Assign simulation parameters from arguments
    tran_step, tran_stop = args.sim_params

    args.ny = args.ny if args.ny is not None else args.nx  # Set NY to NX if not specified
    args.num_vias = args.num_vias if args.num_vias is not None else calculate_default_vias(args.nx, args.ny)
    args.outfile = args.outfile if args.outfile is not None else get_default_outfile(args)

    # Calculate number of digits needed for indices
    digits = len(str(max(args.nx, args.ny)))

    with open(args.outfile, "w") as fout:

        def print_to_file(*args, **kwargs):
            kwargs["file"] = fout
            print(*args, **kwargs)

        # Write header
        print_to_file("* PG sintetica, uniforme \n")

        # Generate resistor network for each layer and capacitors
        total_elements = args.nl * args.nx * args.ny
        with tqdm(total=total_elements, desc="Generating grid elements") as pbar:
            for l in range(1, args.nl + 1):
                for i in range(1, args.nx + 1):
                    for j in range(1, args.ny + 1):
                        # Capacitors
                        print_to_file(
                            f"cL{l}P{i:0{digits}}{j:0{digits}} nL{l}P{i:0{digits}}{j:0{digits}} 0 {capacitor_value}"
                        )
                        if j < args.ny:  # Vertical resistors
                            print_to_file(
                                f"rL{l}P{i:0{digits}}{j:0{digits}}{i:0{digits}}{j+1:0{digits}} nL{l}P{i:0{digits}}{j:0{digits}} nL{l}P{i:0{digits}}{j+1:0{digits}} {resistor_value}"
                            )
                        if i < args.nx:  # Horizontal resistors
                            print_to_file(
                                f"rL{l}P{i:0{digits}}{j:0{digits}}{i+1:0{digits}}{j:0{digits}} nL{l}P{i:0{digits}}{j:0{digits}} nL{l}P{i+1:0{digits}}{j:0{digits}} {resistor_value}"
                            )
                        pbar.update(1)
                    print_to_file()
                print_to_file()

        # Handle VDD on level 1
        print_to_file("* VDD on all nodes on one side on level 1")
        for i in range(1, args.nx + 1):
            print_to_file(
                f"RsrcL1P{i:0{digits}}{1:0{digits}} X_nL1P{i:0{digits}}{1:0{digits}} nL1P{i:0{digits}}{1:0{digits}} {source_resistor_value}"
            )
            print_to_file(f"VL1P{i:0{digits}}{1:0{digits}} X_nL1P{i:0{digits}}{1:0{digits}} 0 {VDD}")
        print_to_file()

        # Add vias between layers if more than one layer
        vias = set()  # Use set for O(1) lookups instead of tuple
        if args.nl > 1:
            print_to_file("* Vias on random nodes on middle of PG")
            with tqdm(total=args.num_vias, desc="Generating vias") as pbar:
                attempts = 0
                max_attempts = args.num_vias * 20  # Prevent infinite loops
                while len(vias) < args.num_vias and attempts < max_attempts:
                    l1 = random.randint(1, args.nl)
                    l2 = l1 - 1 if l1 == args.nl else l1 + 1
                    x = random.randint(2, args.nx)
                    y = random.randint(2, args.ny)
                    via_key = (x, y)  # Simple position key
                    if via_key not in vias:
                        vias.add(via_key)
                        write_via(print_to_file, l1, l2, x, y, via_resistance, digits)
                        pbar.update(1)
                    attempts += 1
        print_to_file()

        # Add current sources
        for i in range(2, args.nx + 1):
            print_to_file(
                f"iB{i:0{digits}}{i:0{digits}} nL1P{i:0{digits}}{i:0{digits}} 0 {current_source_value} pulse({pulse_low}, {pulse_high_value}, {pulse_delay}, {pulse_rise_time}, {pulse_fall_time}, {pulse_on_time}, {pulse_period})"
            )
        print_to_file()

        # Add SPICE control statements
        print_to_file(f".tran {tran_step} {tran_stop}")

        # Print commands for all nodes on level 1 except VDD nodes
        nodes = [f"v(nL1P{i:0{digits}}{j:0{digits}})" for i in range(1, args.nx + 1) for j in range(2, args.ny + 1)]
        print_to_file(".print tran " + " ".join(nodes[-1:]))

        # Write data command
        print_to_file(".end")


if __name__ == "__main__":
    main()
