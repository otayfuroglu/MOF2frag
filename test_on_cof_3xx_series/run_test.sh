#!/bin/bash
# COF-only test runner for UniFrag (fragmentation_oop.py)

PYTHON="/Users/omert/miniconda3/bin/python"
SCRIPT="../fragmentation_oop.py"
COF_DIR="./"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
NC='\033[0m'

RADIUS="${1:-4.0}"

pass=0
fail=0

run_one() {
    local cif="$1"
    local base
    base="$(basename "${cif%.*}")"
    local out="${COF_DIR}/${base}_frag_cof"

    echo -e "${YELLOW}▶ ${base}${NC}"
    output=$($PYTHON "$SCRIPT" "$cif" --kind cof --radius "$RADIUS" --output "$out".xyz 2>&1)
    output=$($PYTHON "$SCRIPT" "$cif" --kind cof --radius "$RADIUS" --minimize --output "$out"_min.xyz 2>&1)

    if [ $? -ne 0 ]; then
        echo -e "  ${RED}✗ FAILED${NC}"
        echo "  $output"
        ((fail++))
        echo ""
        return
    fi

    path_line=$(echo "$output" | grep -E "COF Path A|COF Path B")
    final_line=$(echo "$output" | grep "Final size")

    comp=$($PYTHON -c "from collections import Counter; c=Counter([l.split()[0] for l in open('$out').readlines()[2:] if l.strip()]); print(dict(c))" 2>/dev/null)

    comps=$($PYTHON - <<PY
from collections import deque
import numpy as np

METALS={'Mg','Zn','Cu','Fe','Co','Ni','Mn','Zr','Ti','V','Cr','Al'}
LARGE={'Br','I','S','P','Cl'}

def vb(a,b,d):
    if a in METALS or b in METALS: return d < 2.6
    if 'H' in (a,b): return d < 1.2
    if a in LARGE or b in LARGE: return d < 2.2
    return d < 1.8

sp=[]; co=[]
for ln in open("$out").read().splitlines()[2:]:
    if not ln.strip(): continue
    t=ln.split()
    sp.append(t[0])
    co.append(np.array(list(map(float,t[1:4]))))

n=len(sp)
adj=[[] for _ in range(n)]
for i in range(n):
    for j in range(i+1,n):
        d=np.linalg.norm(co[i]-co[j])
        if vb(sp[i],sp[j],d):
            adj[i].append(j); adj[j].append(i)

vis=set(); sizes=[]
for i in range(n):
    if i in vis: continue
    q=deque([i]); vis.add(i); c=1
    while q:
        u=q.popleft()
        for v in adj[u]:
            if v not in vis:
                vis.add(v); q.append(v); c+=1
    sizes.append(c)

sizes=sorted(sizes, reverse=True)
print(f"{len(sizes)} components; sizes={sizes[:6]}")
PY
)

    echo "  $path_line"
    echo "  $final_line"
    echo "  Composition: $comp"
    echo "  Connectivity: $comps"
    echo -e "  Output: ${GREEN}$out${NC}"
    ((pass++))
    echo ""
}

echo "============================================================"
echo "  UniFrag — COF Test Suite"
echo "============================================================"
echo "  Folder: $COF_DIR"
echo "  Radius: $RADIUS"
echo ""

for cif in "$COF_DIR"/*.cif; do
    run_one "$cif"
done

echo "============================================================"
echo -e "  Results: ${GREEN}$pass passed${NC}, ${RED}$fail failed${NC}"
echo "============================================================"

