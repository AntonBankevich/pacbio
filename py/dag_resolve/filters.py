# from dag_resolve import repeat_graph, sequences
# from dag_resolve.sequences import AlignmentPiece


# def EdgeTransitionFilter(e1, e2, transition_size = 1000):
#     # type: (repeat_graph.Edge, repeat_graph.Edge, int) -> function
#     def f(e1, e2, transition_size, read):
#         # type: (repeat_graph.Edge, repeat_graph.Edge, int, sequences.AlignedRead) -> bool
#         l1 = [] # type: list[AlignmentPiece]
#         l2 = [] # type: list[AlignmentPiece]
#         for al in read.alignments:
#             if al.seg_to.common(e1.suffix(-500)):
#                 l1.append(al)
#             if al.seg_to.common(e2.prefix(500)):
#                 l2.append(al)
#         for al1 in l1:
#             for al2 in l2:
#                 # print al1.__str__(), al2.__str__()
#                 # print al1.rc, al2.rc, al1.seg_from.right - transition_size <= al2.seg_from.left, al2.seg_from.left <= al1.seg_from.right + 50
#                 if al1.rc and al2.rc and al1.seg_from.right - transition_size <= al2.seg_from.left <= al1.seg_from.right + 50:
#                     return True
#                 # print (not al1.rc), (not al2.rc), al2.seg_from.left - transition_size <= al1.seg_from.right, al1.seg_from.right <= al2.seg_from.left + 50
#                 if (not al1.rc) and (not al2.rc) and al2.seg_from.left - transition_size <= al1.seg_from.right <= al2.seg_from.left + 50:
#                     return True
#         return False
#     return lambda read: f(e1, e2, transition_size, read)

