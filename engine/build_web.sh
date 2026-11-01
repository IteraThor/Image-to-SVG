#!/bin/bash
# Build the engine for the web app: web/engine/engine.mjs + engine.wasm (ES module, loaded by
# web/worker.js) and build/engine_node.mjs for command-line checks.
set -euo pipefail
cd "$(dirname "$0")"
# Use an emsdk checked out next to the repo if there is one; otherwise em++ must be on the PATH.
if [ -f ../emsdk/emsdk_env.sh ]; then source ../emsdk/emsdk_env.sh > /dev/null 2>&1; fi
SRC="src/pixel_map.cpp src/pixel_segmentation.cpp src/pixel_segmenter.cpp src/weights_finder.cpp src/subpixel_segmentation.cpp src/subpixel_segmenter.cpp src/shadows.cpp src/subpixel_features.cpp src/boundary_tracer.cpp src/vector_image.cpp src/vector_image_inversion.cpp src/ann_kdtree.cpp src/generative_model.cpp src/optimizer.cpp src/contour_smoother.cpp src/bezier_fitter.cpp src/exporter.cpp src/classifier_trees.cpp src/parameter_classifier.cpp src/palette_finder.cpp src/engine_params.cpp src/engine_params_data.cpp src/pipeline.cpp"
COMMON="-O3 -std=c++20 -Iinclude -sMODULARIZE=1 -sEXPORT_ES6=1 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=67108864 -sSTACK_SIZE=4194304 -sEXPORTED_FUNCTIONS=_isvg_alloc,_isvg_free,_isvg_decode,_isvg_vectorize,_isvg_info,_isvg_export -sEXPORTED_RUNTIME_METHODS=UTF8ToString,stringToUTF8,lengthBytesUTF8,HEAPU8,HEAP32"
mkdir -p ../web/engine build
em++ $COMMON -sENVIRONMENT=web,worker $SRC web/api.cpp -o ../web/engine/engine.mjs
echo "built web/engine/engine.mjs"
em++ $COMMON -sENVIRONMENT=node $SRC web/api.cpp -o build/engine_node.mjs
echo "built build/engine_node.mjs"
