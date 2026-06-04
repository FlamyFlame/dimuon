# A line will be written to args.txt with $pfn replaced by each line in pfns.txt, followed by each set of arguments
# Important: ensure there is a (blank) newline following the last entry in pfns.txt

> args.txt
while read -r pfn; do
  [[ -z "$pfn" ]] && continue
  cat >> args.txt <<EOF
$pfn boo bah blah
$pfn fah fee foe
EOF
done < pfns.txt