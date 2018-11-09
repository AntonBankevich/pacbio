import os
import shutil
import sys

sys.path.append("py")
from common import basic

in_dir = sys.argv[1]
out_dir = sys.argv[2]
basic.ensure_dir_existance(out_dir)
for f in os.listdir(in_dir):
    fig_dir = os.path.join(in_dir, f, "pictures")
    cur = 0
    nums = map(basic.parseNumber, os.listdir(fig_dir))
    last = max(*nums)
    f_name = None
    for f1 in os.listdir(fig_dir):
        if basic.parseNumber(f1) == last:
            f_name = f1
            break
    shutil.copy(os.path.join(fig_dir, f_name), os.path.join(out_dir, f + "-" + f_name))
    
