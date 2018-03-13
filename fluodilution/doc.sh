#!/bin/bash

set -eou pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "doc.sh directory: $DIR"
cd "$DIR"
Rscript scripts/doc.R