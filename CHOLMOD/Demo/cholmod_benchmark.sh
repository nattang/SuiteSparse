cd ..
make demos
cd Demo

#output directory containing A.mtx, B.mtx, X.mtx
OUTPUT_DIR="../../../RXMesh-AMGX/output/"

#if not specified or found, input
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output DIR not found at $OUTPUT_DIR"
    exit 1
fi

for dir in "$OUTPUT_DIR"*/; do
    dir_name=$(basename "$dir")
    dir_name=$(echo "$dir_name" | xargs)
    echo "current directory: $dir"


    A_MATRIX="$dir/A.mtx"
    B_MATRIX="$dir/B.mtx"
    X_MATRIX="$dir/X.mtx"

    if [ ! -f "$A_MATRIX" ] || [ ! -f "$B_MATRIX" ] || [ ! -f "$X_MATRIX" ]; then
        echo "Skipping $dir_name - Doesn't contain necessary .mtx files"
        continue
    fi

    echo "Running CHOLMOD for $dir_name..."
        CHOLMOD_USE_GPU=1 ../build/cholmod_dl_simple \
        -A $A_MATRIX \
        -B $B_MATRIX \
        -EX $X_MATRIX
    echo "-------------------------------"
done