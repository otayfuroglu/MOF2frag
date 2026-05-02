#!/bin/bash
# Quick regression test for MOF2frag (one MOF per family)

PYTHON="/Users/omert/miniconda3/bin/python"
OOP_SCRIPT="$(dirname "$0")/fragmentation_oop.py"
ENGINE="oop"
KIND="mof"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
NC='\033[0m'

pass=0
fail=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --engine)
            ENGINE="$2"
            shift 2
            ;;
        --kind)
            KIND="$2"
            shift 2
            ;;
        *)
            echo "Unknown arg: $1"
            echo "Usage: $0 [--engine oop] [--kind mof|cof]"
            exit 1
            ;;
    esac
done

if [[ "$ENGINE" != "oop" ]]; then
    echo "Invalid --engine: $ENGINE (only 'oop' is supported)"
    exit 1
fi

run_test() {
    local label="$1"
    local cif="$2"
    local radius="$3"
    local nmetals="$4"
    local extra_args="$5"

    local suffix="_frag_fast"
    [[ "$extra_args" == *"--minimize"* ]] && suffix="_frag_fast_min"
    local outfile="$(dirname "$cif")/$(basename "${cif%.*}")${suffix}.xyz"

    echo -e "${YELLOW}▶ $label${NC}"
    if [[ "$ENGINE" == "oop" ]]; then
        output=$($PYTHON "$OOP_SCRIPT" "$cif" --kind "$KIND" --radius "$radius" --nmetals "$nmetals" --output "$outfile" $extra_args 2>&1)
    else
        echo "  [skip] legacy engine removed; use --engine oop"
        return
    fi

    if [ $? -ne 0 ]; then
        echo -e "  ${RED}✗ FAILED${NC}"
        echo "  $output"
        ((fail++))
        echo ""
        return
    fi

    path_line=$(echo "$output" | grep -E "Path A|Path B|Path C")
    final_line=$(echo "$output" | grep "Final size")
    echo "  $path_line"
    echo "  $final_line"
    echo -e "  Output: ${GREEN}$outfile${NC}"
    ((pass++))
    echo ""
}

echo "============================================================"
echo "  MOF2frag — Fast Test Suite (1 per family)"
echo "============================================================"
echo "  Engine: $ENGINE"
[[ "$ENGINE" == "oop" ]] && echo "  Kind: $KIND"
echo ""

# One representative MOF per family
MOF74_CIF="$(dirname "$0")/test/1516648.cif"
IRMOF_CIF="$(dirname "$0")/test_on_irmof_series/IRMOF-1.cif"
MIL_CIF="$(dirname "$0")/test_on_mil_series/MIL-100.cif"
PCN_CIF="$(dirname "$0")/test_on_pcn_series/PCN-60.cif"

echo -e "${CYAN}Normal mode${NC}"
run_test "MOF-74: 1516648" "$MOF74_CIF" 4.0 3 ""
run_test "IRMOF: IRMOF-1" "$IRMOF_CIF" 4.0 3 ""
run_test "MIL: MIL-100" "$MIL_CIF" 4.0 3 ""
run_test "PCN: PCN-60" "$PCN_CIF" 4.0 3 ""

echo -e "${CYAN}Minimize mode${NC}"
run_test "MOF-74 [min]: 1516648" "$MOF74_CIF" 4.0 3 "--minimize"
run_test "IRMOF [min]: IRMOF-1" "$IRMOF_CIF" 4.0 3 "--minimize"
run_test "MIL [min]: MIL-100" "$MIL_CIF" 4.0 3 "--minimize"
run_test "PCN [min]: PCN-60" "$PCN_CIF" 4.0 3 "--minimize"

echo "============================================================"
echo -e "  Results: ${GREEN}$pass passed${NC}, ${RED}$fail failed${NC}"
echo "============================================================"
