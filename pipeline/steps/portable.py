SI="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/supplementary/build_supplementary_tables.py"
s=open(SI,encoding="utf-8").read()
# 插入可移植 BASE 检测（仅插一次）
if "BASE_DATA" not in s:
    inject=('import os\n'
            'BASE_DATA = next((b for b in ("/mnt/data2/墨江紫米研究","/Volumes/data2/墨江紫米研究") if os.path.exists(b)), "/mnt/data2/墨江紫米研究")\n')
    s=s.replace("import os\n", inject, 1)
# 把绝对前缀替换为 BASE_DATA 拼接
s=s.replace('"/Volumes/data2/墨江紫米研究', 'BASE_DATA + "')
open(SI,"w",encoding="utf-8").write(s)
print("done; /Volumes refs left:", s.count("/Volumes/data2"))
