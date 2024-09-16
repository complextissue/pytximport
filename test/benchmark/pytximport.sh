rm *.json
python3 bench.py -o ./pytximport_time.json -p 1 -n 11 -w 0
python3 -m pyperf stats ./pytximport_time.json
python3 bench.py -o ./pytximport_memory.json -p 1 -n 11 -w 0 --track-memory
python3 -m pyperf stats ./pytximport_memory.json
