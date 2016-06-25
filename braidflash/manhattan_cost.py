# Pranav Gokhale
# Usage: python arrange.py input_trace_file.tr [width height]
# (if width and height are not specified, default to sidelength
#  of smallest square big enough to contain all nodes)
# Input file: .tr file
# Output:
#    - placement of nodes on the rectangle

import sys, math

def cost_for_trace(tracefile):
  """Converts .tr file to .graph file."""
  assert '.tr' in tracefile
  tracefile = open(tracefile)
  total = 0
  module_cost = -1

  for line in tracefile:
    if line.startswith('module:'):
      if module_cost >= 0:
        total += module_cost
        print module_cost
      module_cost = 0
      tokens = line.split(' ')
      module_name = tokens[1]
      print "Total Manhattan Cost for ", module_name,
      continue
    elif line.startswith('num_nodes:'):
      tokens = line.split(' ')
      num_nodes = int(tokens[1])
      square_size = _square_size(num_nodes)
      continue
    elif not line.startswith('ID:'):  # skip the header lines
      continue
    elif 'DST:' not in line:          # skip single qubit gates
      continue

    tokens = line.split(' ')
    assert tokens[4].startswith('SRC:')
    assert tokens[6].startswith('DST:')
    src = int(tokens[5])
    dst = int(tokens[7])
    assert src < num_nodes
    assert dst < num_nodes

    edge_cost_y = abs(src / square_size - dst / square_size)
    edge_cost_x = abs(src % square_size - dst % square_size)
    module_cost += (edge_cost_x + edge_cost_y)

  total += module_cost
  print module_cost  # last one won't get printed by loop  
  print '\n\nTotal cost for all modules (i.e. sum of costs above) is:', total

def _square_size(module_num_nodes):
  return int(math.ceil(math.sqrt(module_num_nodes)))


if __name__ == "__main__":
  cost_for_trace(sys.argv[1])
