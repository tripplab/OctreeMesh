#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print("Usage: generate_hex_grid.py <n_hex> <output.post.msh>")
    sys.exit(2)

n_hex = int(sys.argv[1])
out = sys.argv[2]

with open(out, "w", encoding="utf-8") as f:
    f.write('MESH "hex" dimension 3 ElemType Hexahedra Nnode 8\n\n')
    f.write('Coordinates\n')
    node_id = 1
    elem_id = 1
    for e in range(n_hex):
      x = float(e)
      coords = [
        (x,0,0),(x+1,0,0),(x+1,1,0),(x,1,0),
        (x,0,1),(x+1,0,1),(x+1,1,1),(x,1,1)
      ]
      for c in coords:
        f.write(f"{node_id} {c[0]} {c[1]} {c[2]}\n")
        node_id += 1
    f.write('End Coordinates\n')
    f.write('Elements\n')
    nid = 1
    for e in range(n_hex):
      f.write(f"{elem_id} {nid} {nid+1} {nid+2} {nid+3} {nid+4} {nid+5} {nid+6} {nid+7}\n")
      elem_id += 1
      nid += 8
    f.write('End Elements\n')
