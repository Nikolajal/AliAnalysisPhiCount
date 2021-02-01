# Execute alidock within the appropriate Python virtual environment
function alidock() {( source "/home/nrubini/.virtualenvs/alidock/bin/activate" && command alidock "$@"; exit $?; )}
