import os
import shutil
from tqdm import tqdm

def move_dumpall_files_2layers(source_root, target_root):
    dumpall_dirs = []
    max_depth = 1
    root_depth = source_root.rstrip(os.sep).count(os.sep)

    # 2層目までのdumpall.datを探す
    for dirpath, dirnames, filenames in os.walk(source_root):
        current_depth = dirpath.count(os.sep) - root_depth
        if current_depth > max_depth:
            # ディレクトリの走査を止める
            dirnames[:] = []  # 現在のdirpath以下を探索しない
            continue
        if 'dumpall.dat' in filenames:
            dumpall_dirs.append(dirpath)

    # tqdmで進捗表示
    for dirpath in tqdm(dumpall_dirs, desc="Moving dumpall.dat files", unit="file"):
        relative_path = os.path.relpath(dirpath, source_root)
        target_path = os.path.join(target_root, relative_path)
        os.makedirs(target_path, exist_ok=True)

        source_file = os.path.join(dirpath, 'dumpall.dat')
        target_file = os.path.join(target_path, 'dumpall.dat')
        shutil.move(source_file, target_file)

# 使用例
source_root = "H:/hata/test_compare"
target_root = "H:/hata2025/1332_120_100"
move_dumpall_files_2layers(source_root, target_root)
