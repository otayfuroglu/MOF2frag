#!/bin/bash
# PCN-only test runner for fragmentation_oop.py

PYTHON="/Users/omert/miniconda3/bin/python"
SCRIPT="../fragmentation_oop.py"
PCN_DIR="./"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
NC='\033[0m'

pass=0
fail=0

run_one() {
    local cif="$1"
    local mode="$2"      # norm|min
    local extra_args="$3"

    local base
    base="$(basename "${cif%.*}")"

    local outfile="${PCN_DIR}/${base}_frag_oop_pcn.xyz"
    [[ "$mode" == "min" ]] && outfile="${PCN_DIR}/${base}_frag_oop_pcn_min.xyz"

    echo -e "${YELLOW}▶ ${base} (${mode})${NC}"
    # output=$($PYTHON "$SCRIPT" "$cif" --kind mof --radius 4.0 --nmetals 3 --output "$outfile" $extra_args 2>&1)
    output=$($PYTHON "$SCRIPT" "$cif" --radius 4.0 --nmetals 3 --output "$outfile" $extra_args 2>&1)

    if [ $? -ne 0 ]; then
        echo -e "  ${RED}✗ FAILED${NC}"
        echo "  $output"
        ((fail++))
        echo ""
        return
    fi

    path_line=$(echo "$output" | grep -E "Path A|Path B|Path C")
    final_line=$(echo "$output" | grep "Final size")
    comp=$($PYTHON -c "from collections import Counter; c=Counter([l.split()[0] for l in open('$outfile').readlines()[2:] if l.strip()]); print(dict(c))" 2>/dev/null)

    echo "  $path_line"
    echo "  $final_line"
    echo "  Composition: $comp"
    echo -e "  Output: ${GREEN}$outfile${NC}"
    ((pass++))
    echo ""
}

echo "============================================================"
echo "  PCN-only OOP Test Suite (fragmentation_oop.py)"
echo "============================================================"
echo ""

for cif in "$PCN_DIR"/*.cif; do
    run_one "$cif" "norm" ""
done

echo -e "${CYAN}Minimize mode${NC}"
for cif in "$PCN_DIR"/*.cif; do
    run_one "$cif" "min" "--minimize"
done

echo "============================================================"
echo -e "  Results: ${GREEN}$pass passed${NC}, ${RED}$fail failed${NC}"
echo "============================================================"
