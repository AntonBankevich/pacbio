import os
import sys


def choose_datasets(dir):
    list_file = os.path.join(dir, "list.txt")
    list = [s.strip() for s in open(list_file, "r").readlines()]
    for dataset_id in list:
        assembly_result = os.path.join(dir, "assembly", dataset_id, "lines.info")
        if not os.path.isfile(assembly_result):
            print dataset_id
            continue
        ok = True
        for s in open(assembly_result).readlines():
            s = s.strip()
            if s.startswith("Line:") and (s.find("(") == -1 or s.find(")") == -1):
                ok = False
                break
        if ok:
            print dataset_id

if __name__ == "__main__":
    choose_datasets(sys.argv[1])