# BLAST_DBS_DOWNLOADER
A threaded python script that automates download and extraction of NCBI Blast databases

Execution:  python update_blastdbs.py <database names> <threads>, interactively guides through do or update

example:  python update_blastdbs.py nt nr -j 4

If you execute this script with no parameters it will bring down a complete list of available databases.
