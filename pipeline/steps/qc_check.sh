M="/mnt/data2/墨江紫米研究/12_论文Paper1_Manuscript/manuscript.tex"
echo "=== cite keys vs bibitem keys ==="
grep -oE '\\cite\{[^}]*\}' "$M" | sed 's/\\cite{//;s/}//' | tr ',' '\n' | sort -u > /tmp/cited.txt
grep -oE '\\bibitem\{[^}]*\}' "$M" | sed 's/\\bibitem{//;s/}//' | sort -u > /tmp/defined.txt
echo "cited 但未定义 (应空):"; comm -23 /tmp/cited.txt /tmp/defined.txt
echo "定义但未引用:"; comm -13 /tmp/cited.txt /tmp/defined.txt
echo "cite 总数 / bibitem 总数: $(wc -l < /tmp/cited.txt) / $(wc -l < /tmp/defined.txt)"
echo ""
echo "=== 'copy-number-expanded' 上下文 ==="
grep -n "copy-number-expanded" "$M"
echo ""
echo "=== 'will be deposited' 上下文 (Data availability) ==="
grep -n -A2 "will be deposited" "$M"
echo ""
echo "=== 哪些 Table S 在正文被引用 (S1-S11) ==="
for n in 1 2 3 4 5 6 7 8 9 10 11; do printf "S%s:" $n; grep -oE "Table~?S$n([^0-9]|\$)" "$M" | head -1 | grep -q . && echo " cited" || echo " NOT cited"; done
