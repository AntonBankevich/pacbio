import sys


sys.path.append("py")
from common import basic, SeqIO
from common.SimpleGraph import SimpleGraph

graph = SimpleGraph().ReadGFA(sys.argv[1])
for e_id in graph.e:
    if basic.isCanonocal(e_id):
        SeqIO.write(graph.e[e_id], sys.stdout, "fasta")
