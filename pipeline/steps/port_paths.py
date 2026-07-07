import os,glob,re
V="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/V2"
old_data="/sessions/friendly-exciting-dirac/mnt/墨江紫米研究/"
new_data="/mnt/data2/墨江紫米研究/"
old_out="/sessions/friendly-exciting-dirac/mnt/outputs"
new_out=V
files=glob.glob(V+"/*.py")+glob.glob(V+"/figs_v2/*.py")
n=0
for f in files:
    s=open(f,encoding="utf-8").read(); o=s
    s=s.replace(old_data,new_data).replace(old_out,new_out)
    if s!=o:
        open(f+".orig","w",encoding="utf-8").write(o)  # backup
        open(f,"w",encoding="utf-8").write(s); n+=1
        print("ported:",os.path.basename(f))
print(f"{n} files ported (.orig backups kept)")
# verify no stale session paths remain
import subprocess
left=subprocess.run(["grep","-rl","friendly-exciting-dirac",V],capture_output=True,text=True).stdout.strip()
print("remaining session-path files:", left.count("\n")+1 if left else 0)
