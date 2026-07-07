SI="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/supplementary/build_supplementary_tables.py"
s=open(SI,encoding="utf-8").read()
bad='BASE_DATA = next((b for b in ("/mnt/data2/墨江紫米研究",BASE_DATA + "") if os.path.exists(b)), "/mnt/data2/墨江紫米研究")'
good='BASE_DATA = next((b for b in ("/mnt/data2/墨江紫米研究","/Volumes/data2/墨江紫米研究") if os.path.exists(b)), "/mnt/data2/墨江紫米研究")'
assert s.count(bad)==1, "bad-line count "+str(s.count(bad))
open(SI,"w",encoding="utf-8").write(s.replace(bad,good))
print("fixed BASE_DATA line")
