import os
from Bio import SeqIO
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm  # 🚀 加载进度条

# ==== 参数设置 ====
genome_file = "Anoxybacillus_gonensis.fasta"
short_rna_file = "HYER.fasta"
output_file = "foldalign_all_windows_output.txt"
window_size = 700
step_size = 50
num_workers = 4  # 改成你电脑的核心数

temp_dir = "temp_windows"
os.makedirs(temp_dir, exist_ok=True)

# ==== 读取长基因组 ====
genome_record = next(SeqIO.parse(genome_file, "fasta"))
genome_seq = str(genome_record.seq)

# ==== 函数：处理一个窗口 ====
def process_window(i):
    window_seq = genome_seq[i:i + window_size]
    window_filename = os.path.join(temp_dir, f"window_{i}.fasta")

    with open(window_filename, "w") as w_f:
        w_f.write(f">window_{i}_{i+window_size}\n{window_seq}\n")

    cmd = [
        "/Users/yss/Desktop/1_R/1_101_Bio_Info/foldalign.2.5.3/bin/foldalign",
        "-max_length", "700",
        "-output_format", "summary",
        "-number_of_processors", "1",
        short_rna_file,
        window_filename
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)  # 加个timeout保险
        if result.returncode != 0:
            return i, f"[Foldalign 运行失败] stderr:\n{result.stderr}"
        return i, result.stdout
    except Exception as e:
        return i, f"[异常错误] 窗口{i}: {str(e)}"

if __name__ == '__main__':
    import shutil

# ==== 并行处理窗口，比对过程可视化 ====
    all_positions = list(range(0, len(genome_seq) - window_size + 1, step_size))
    results = []

    print(f"开始运行 Foldalign 并行比对（共 {len(all_positions)} 个窗口）...")

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_window, i): i for i in all_positions}
        for future in tqdm(as_completed(futures), total=len(futures), desc="比对进度"):
            i = futures[future]
            try:
                index, output = future.result()
                results.append((index, output))
            except Exception as e:
                #print(f"[错误] 窗口 {i}: {e}")
                a = 1

    # ==== 保存输出 ====
    results.sort(key=lambda x: x[0])  # 按窗口顺序排列

    with open(output_file, "w") as out_f:
        for index, output in results:
            out_f.write(f"\n>>> Window {index}-{index + window_size} <<<\n")
            out_f.write(output)

    print("✅ 所有窗口比对完成，结果保存在：", output_file)

    # ==== 清理临时文件 ====
    shutil.rmtree(temp_dir)
