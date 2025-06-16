#!/usr/bin/env python3
"""Power grid circuit generator."""

import sys
import random
import argparse
from pathlib import Path
import math

# Voltage values
VDD = 1.8
GND = 0.0

def write_via(out_print, l1, l2, x, y, val, digits):
    """Write a via (vertical connection) between two layers."""
    out_print(f"VL{l1}P{x:0{digits}}{y:0{digits}}L{l2}P{x:0{digits}}{y:0{digits}} nL{l1}P{x:0{digits}}{y:0{digits}} nL{l2}P{x:0{digits}}{y:0{digits}} {val}")

def calculate_default_vias(nx, ny):
    """Calculate default number of vias as 20% of grid size."""
    return int(0.2 * nx * ny)

def get_default_outfile(args):
    """Generate default output filename based on grid parameters."""
    return f"grid_{args.nx}x{args.ny}_l{args.nl}_v{args.num_vias}.sp"

def main():
    parser = argparse.ArgumentParser(description='Generate a power grid netlist.')
    parser.add_argument('outfile', nargs='?', help='Output file path (default: grid_NXxNY_lNL_vNVIAS.sp)')
    parser.add_argument('--nx', type=int, default=4, help='Grid size in X dimension (default: 4)')
    parser.add_argument('--ny', type=int, default=None, help='Grid size in Y dimension (default: nx)')
    parser.add_argument('--nl', type=int, default=2, help='Number of layers (default: 2)')
    parser.add_argument('--num-vias', type=int, help='Number of vias to add (default: 20%% of grid size)')
    
    args = parser.parse_args()
    args.ny = args.ny if args.ny is not None else args.nx  # Set NY to NX if not specified
    args.num_vias = args.num_vias if args.num_vias is not None else calculate_default_vias(args.nx, args.ny)
    args.outfile = args.outfile if args.outfile is not None else get_default_outfile(args)
    
    # Calculate number of digits needed for indices
    digits = len(str(max(args.nx, args.ny)))
    
    with open(args.outfile, "w") as fout:
        def print_to_file(*args, **kwargs):
            kwargs['file'] = fout
            print(*args, **kwargs)
        
        # Write header
        print_to_file("* PG sintetica, uniforme \n")
        
        # Generate resistor network for each layer and capacitors
        for l in range(1, args.nl + 1):
            for i in range(1, args.nx + 1):
                for j in range(1, args.ny + 1):
                    # Capacitors
                    print_to_file(f"cL{l}P{i:0{digits}}{j:0{digits}} nL{l}P{i:0{digits}}{j:0{digits}} 0 5.37096111111111e-11")
                    if j < args.ny:  # Vertical resistors
                        print_to_file(f"rL{l}P{i:0{digits}}{j:0{digits}}{i:0{digits}}{j+1:0{digits}} nL{l}P{i:0{digits}}{j:0{digits}} nL{l}P{i:0{digits}}{j+1:0{digits}} 2.50e-01")
                    if i < args.nx:  # Horizontal resistors
                        print_to_file(f"rL{l}P{i:0{digits}}{j:0{digits}}{i+1:0{digits}}{j:0{digits}} nL{l}P{i:0{digits}}{j:0{digits}} nL{l}P{i+1:0{digits}}{j:0{digits}} 2.50e-01")
                print_to_file()
            print_to_file()
        
        # Handle VDD on level 1
        print_to_file("* VDD on all nodes on one side on level 1")
        for i in range(1, args.nx + 1):
            print_to_file(f"RsrcL1P{i:0{digits}}{1:0{digits}} X_nL1P{i:0{digits}}{1:0{digits}} nL1P{i:0{digits}}{1:0{digits}} 2.5e-02")
            print_to_file(f"VL1P{i:0{digits}}{1:0{digits}} X_nL1P{i:0{digits}}{1:0{digits}} 0 {VDD}")
        print_to_file()
        
        # Add vias between layers if more than one layer
        vias = ()
        if args.nl > 1:
            print_to_file("* Vias on random nodes on middle of PG")
            for _ in range(args.num_vias):
                l1 = random.randint(1, args.nl)
                l2 = l1 - 1 if l1 == args.nl else l1 + 1
                x = min(random.randint(2, args.nx + 1), args.nx)
                y = min(random.randint(2, args.ny + 1), args.ny)
                if f"{x:0{digits}}/{y:0{digits}}" not in vias:
                    write_via(print_to_file, l1, l2, x, y, GND, digits)
                    vias += (f"{x:0{digits}}/{y:0{digits}}", )
        print_to_file()
        
        # Add current sources
        for i in range(2, args.nx + 1):
            print_to_file(f"iB{i:0{digits}}{i:0{digits}} nL1P{i:0{digits}}{i:0{digits}} 0 2.5099e-05 pulse(2.5099e-05, 0.0627475, 1.0e-9, 1e-19, 1e-10, 1e-19, 3e-09)")
        print_to_file()

        # Add SPICE control statements
        print_to_file(".tran 1.000000001e-12 1e-8")
        
        # Print commands for all nodes on level 1 except VDD nodes
        nodes = [f"v(nL1P{i:0{digits}}{j:0{digits}})" for i in range(1, args.nx + 1) for j in range(2, args.ny + 1)]
        print_to_file(".print tran " + " ".join(nodes))
        
        # Write data command
        print_to_file(".end")

if __name__ == "__main__":
    main()
