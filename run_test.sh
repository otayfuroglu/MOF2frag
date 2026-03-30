#!/bin/bash
# ============================================================
# MOF2frag Test Script
# Tests Path A (discrete 0D clusters) and Path B (infinite rods)
# ============================================================

PYTHON="/Users/omert/miniconda3/bin/python"
SCRIPT="$(dirname "$0")/fragmenter.py"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
NC='\033[0m'

echo "============================================================"
echo "  MOF2frag — Automated Test Suite"
echo "============================================================"
echo ""

pass=0
fail=0

run_test() {
    local label="$1"
    local cif="$2"
    local radius="$3"
    local nmetals="$4"
    local extra_args="$5"
    local suffix="_frag"
    [[ "$extra_args" == *"--minimize"* ]] && suffix="_frag_min"
    local outfile="$(dirname "$cif")/$(basename "${cif%.*}")${suffix}.xyz"

    echo -e "${YELLOW}▶ $label${NC}"
    output=$($PYTHON "$SCRIPT" "$cif" --radius "$radius" --nmetals "$nmetals" --output "$outfile" $extra_args 2>&1)

    if [ $? -ne 0 ]; then
        echo -e "  ${RED}✗ FAILED to run${NC}"
        echo "  $output"
        ((fail++))
        echo ""
        return
    fi

    sbu_line=$(echo "$output" | grep -E "Path A|Path B|Path C")
    final_line=$(echo "$output" | grep "Final size")
    prune_line=$(echo "$output" | grep "Removed")
    comp=$($PYTHON -c "from collections import Counter; c=Counter([l.split()[0] for l in open('$outfile').readlines()[2:] if l.strip()]); print(dict(c))" 2>/dev/null)

    echo "  $sbu_line"
    [ -n "$prune_line" ] && echo "  $prune_line"
    echo "  $final_line"
    echo "  Composition: $comp"
    echo -e "  Output: ${GREEN}$outfile${NC}"
    ((pass++))
    echo ""
}

# ===== PATH B: MOF-74 series (infinite metal rods) =====
echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  PATH B — Infinite Metal Lattice (MOF-74 series)${NC}"
echo -e "${CYAN}============================================================${NC}"
echo ""

for cif in "$(dirname "$0")"/test/*.cif; do
    name=$(basename "$cif" .cif)
    run_test "MOF-74: $name (3 metals)" "$cif" 4.0 3 ""
done

# Also test with 5 metals on Cu MOF-74
if [ -f "$(dirname "$0")/test/1516648.cif" ]; then
    run_test "MOF-74: 1516648 (5 metals)" "$(dirname "$0")/test/1516648.cif" 4.0 5 ""
fi

# ===== PATH A: IRMOF series (discrete 0D Zn4O clusters) =====
echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  PATH A — Discrete 0D Clusters (IRMOF series)${NC}"
echo -e "${CYAN}============================================================${NC}"
echo ""

for cif in "$(dirname "$0")"/test_on_irmof_series/*.cif; do
    name=$(basename "$cif" .cif)
    run_test "IRMOF: $name" "$cif" 4.0 3 ""
done

# ===== PATH A & B: MIL Series (Discrete & Infinite chains) =====
echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  PATH A & PATH B — MIL Series (Discrete & Infinite)${NC}"
echo -e "${CYAN}============================================================${NC}"
echo ""

for cif in "$(dirname "$0")"/test_on_mil_series/*.cif; do
    name=$(basename "$cif" .cif)
    run_test "MIL: $name" "$cif" 4.0 3 ""
done

# ===== PATH A: PCN Series (Discrete clusters) =====
echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  PATH A — PCN Series${NC}"
echo -e "${CYAN}============================================================${NC}"
echo ""

for cif in "$(dirname "$0")"/test_on_pcn_series/*.cif; do
    name=$(basename "$cif" .cif)
    run_test "PCN: $name" "$cif" 4.0 3 ""
done

# ============================================================
#  MINIMIZE MODE — All series re-run with --minimize
# ============================================================

echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  MINIMIZE MODE — MOF-74 Series${NC}"
echo -e "${CYAN}============================================================${NC}"
echo ""

for cif in "$(dirname "$0")"/test/*.cif; do
    name=$(basename "$cif" .cif)
    run_test "MOF-74 [min]: $name (3 metals)" "$cif" 4.0 3 "--minimize"
done

echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  MINIMIZE MODE — IRMOF Series${NC}"
echo -e "${CYAN}============================================================${NC}"
echo ""

for cif in "$(dirname "$0")"/test_on_irmof_series/*.cif; do
    name=$(basename "$cif" .cif)
    run_test "IRMOF [min]: $name" "$cif" 4.0 3 "--minimize"
done

echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  MINIMIZE MODE — MIL Series${NC}"
echo -e "${CYAN}============================================================${NC}"
echo ""

for cif in "$(dirname "$0")"/test_on_mil_series/*.cif; do
    name=$(basename "$cif" .cif)
    run_test "MIL [min]: $name" "$cif" 4.0 3 "--minimize"
done

echo -e "${CYAN}============================================================${NC}"
echo -e "${CYAN}  MINIMIZE MODE — PCN Series${NC}"
echo -e "${CYAN}============================================================${NC}"
echo ""

for cif in "$(dirname "$0")"/test_on_pcn_series/*.cif; do
    name=$(basename "$cif" .cif)
    run_test "PCN [min]: $name" "$cif" 4.0 3 "--minimize"
done

# ===== Summary =====
echo "============================================================"
echo -e "  Results: ${GREEN}$pass passed${NC}, ${RED}$fail failed${NC}"
echo "  Output XYZ files saved next to their source CIF files."
echo "============================================================"
