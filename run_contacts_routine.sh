#!/bin/bash  
# script to create all minimal distances with support for multiple trajectories

# Print usage instructions
usage() {
    echo "Usage: $0 -traj <path> -pdb <file> -min <val> -max <val> -ndx <file> -sys <name> -traj_mode <mode> -mode <mode>"
    echo "  -traj: MDAnalysis-supported trajectory file OR trajectory folder"
    echo "  -pdb: pdb file"
    echo "  -min: minimum threshold (0-1)"
    echo "  -max: maximum threshold (0-1)"
    echo "  -ndx: indices file"
    echo "  -sys: base name for system, for output files"
    echo "  -traj_mode: traj mode [multi|single] trajectories"
    echo "  -mode: threshold mode [overall|per-trajectory]"
    echo "           overall: contacts averaged across all trajectories"
    echo "           per-trajectory: contacts must meet threshold in each trajectory"
    echo ""
    echo "Examples:"
    echo "  Single trajectory:"
    echo "    $0 -traj traj.dcd -pdb system.pdb -min 0.3 -max 0.9 -ndx indices.ndx -sys system_name -traj_mode single -mode overall"
    echo ""
    echo "  Multiple trajectories (overall mode):"
    echo "    $0 -traj /path/to/traj_folder -pdb system.pdb -min 0.3 -max 0.9 -ndx indices.ndx -sys system_name -traj_mode multi -mode overall"
    echo ""
    echo "  Multiple trajectories (per-trajectory mode):"
    echo "    $0 -traj /path/to/traj_folder -pdb system.pdb -min 0.3 -max 0.9 -ndx indices.ndx -sys system_name -traj_mode multi -mode per-trajectory"
    echo ""
    echo "  A folder may contain any mixture of MDAnalysis-supported trajectory formats."
} 

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -traj|-xtc) TRAJECTORY="$2"; shift ;;
        -pdb) PDB="$2"; shift ;;
        -min) MIN_THR="$2"; shift ;;
        -max) MAX_THR="$2"; shift ;;
        -ndx) INDEX="$2"; shift ;;
        -sys) SYSTEM="$2"; shift ;;
        -traj_mode) TRAJ_MODE="$2"; shift ;;
        -mode) MODE="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
    shift
done

# Check if all required arguments are provided
if [ -z "$TRAJECTORY" ] || [ -z "$PDB" ] || [ -z "$MIN_THR" ] || [ -z "$MAX_THR" ] || [ -z "$INDEX" ] || [ -z "$SYSTEM" ] || [ -z "$TRAJ_MODE" ] || [ -z "$MODE" ]; then
    echo "Error: Missing required arguments."
    usage
    exit 1
fi

# Validate trajectory mode
if [ "$TRAJ_MODE" != "single" ] && [ "$TRAJ_MODE" != "multi" ]; then
    echo "Error: -traj_mode must be 'single' or 'multi'"
    usage
    exit 1
fi

# Validate contact mode
if [ "$MODE" != "overall" ] && [ "$MODE" != "per-trajectory" ]; then
    echo "Error: -mode must be 'overall' or 'per-trajectory'"
    usage
    exit 1
fi  

if [ "$TRAJ_MODE" == "single" ] && [ "$MODE" == "per-trajectory" ]; then
    echo "Error: threshold mode 'per-trajectory' is only allowed with trajectory mode 'multi'"
    usage
    exit 1
fi

# Check if multi-trajectory mode
if [ "$TRAJ_MODE" == "multi" ]; then
    TRAJ_ARG="--trajectory-list"
    echo "Running in multi-trajectory mode"
    echo "Trajectory folder: $TRAJECTORY"
    echo "Mode: $MODE"
    if [ "$MODE" == "overall" ]; then
        echo "threshold is set over the total length of simulated data"
    elif [ "$MODE" == "per-trajectory" ]; then
        echo "threshold is set over each single trajectory"
    fi
else
    TRAJ_ARG="--trajectory"
    echo "Running in single-trajectory mode"
    echo "Trajectory file: $TRAJECTORY"
fi  

if [ ! -f "$PDB" ] || [ ! -f "$INDEX" ]; then
    echo "Error: PDB or INDEX file does not exist"
    usage
    exit 1
fi

# Check the trajectory path based on mode
if [ "$TRAJ_MODE" == "multi" ]; then
    if [ ! -d "$TRAJECTORY" ]; then
        echo "Error: In multi mode, -traj must be a trajectory directory"
        usage
        exit 1
    fi
else
    if [ ! -f "$TRAJECTORY" ]; then
        echo "Error: In single mode, -traj must be a trajectory file"
    usage
    exit 1
    fi
fi 

# created files  
THR_SUFFIX="${MIN_THR}-${MAX_THR}"  
IS_MINDIST="${SYSTEM}.is_mindist"
ATOMDIST="${SYSTEM}.all_thr${THR_SUFFIX}_selected_atom_distances"

# estimate all mindist
if [[ ! -e "${IS_MINDIST}" ]]; then
    echo ""
    echo "step 1 - estimate_contacts.py: compute contacts"
    time python estimate_contacts.py \
        --top $PDB \
        "$TRAJ_ARG" "$TRAJECTORY" \
        --index $INDEX \
        --output $IS_MINDIST \
        --mode $MODE \
        --verbose
else
    echo "skipping step 1 (${IS_MINDIST} already exists)"
fi

# select formed mindist
echo ""
echo "step 2 - extract_indices.py: select residue pairs where mindist is contact for more than min threshold time and less than max threshold time"
time python extract_indices.py \
    --is-contacts $IS_MINDIST \
    --threshold $MIN_THR \
    --max-threshold $MAX_THR 

# extract all atom pairwise distances of selected residues
echo ""
echo "step 3 - contacts.py: extract all atom pairwise dists between residues in .ndx file"
time python contacts.py \
    "$TRAJ_ARG" "$TRAJECTORY" \
    -s $PDB \
    -n ${IS_MINDIST}.thr${THR_SUFFIX}.ndx \
    -o $ATOMDIST

# extract minimal distances between all atom pairs forming a contact more often
# than the the given threshold
echo ""
echo "step 4 - extract_contacts.py: extracting final minimal distances"
time python extract_contacts.py \
    --contacts $ATOMDIST \
    --index ${ATOMDIST}.atom_indices \
    --threshold ${MIN_THR} \
    --output ${SYSTEM}.mindist  

echo ""
echo "Final output files:"
echo "  - ${SYSTEM}.mindist (final minimal distances)"
echo "  - ${SYSTEM}.mindist.ndx (indices of selected residue pairs)"
#end
