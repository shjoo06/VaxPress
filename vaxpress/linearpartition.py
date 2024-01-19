import os
import sys
import subprocess as sp
from tqdm import tqdm
from select import select
from fcntl import fcntl, F_GETFL, F_SETFL

def read_live_updates(stdout, bufsize=8192):
    fd = stdout.fileno()
    flags = fcntl(fd, F_GETFL)
    fcntl(fd, F_SETFL, flags | os.O_NONBLOCK)

    buf = []

    while True:
        r, _, _ = select([fd], [], [], 0.1)
        if not r:
            continue

        data = os.read(fd, bufsize)
        if not data:
            break

        if b'\r' in data:
            chunks = data.split(b'\r')
            yield (b''.join(buf) + chunks[0]).decode()

            for chunk in chunks[1:-1]:
                yield chunk.decode()

            buf[:] = [chunks[-1]]
        else:
            buf.append(data)

    if buf:
        yield b''.join(buf).decode()

def run_linearpartition(linearpartition_dir, sequence):
    linearpartition_bin = os.path.join(linearpartition_dir, 'linearpartition')
    linearpartition_bin = os.path.abspath(linearpartition_bin)
    if not os.path.exists(linearpartition_bin):
        raise FileNotFoundError('linearpartition binary not found at {}'.format(linearpartition_bin))

    # store LP output file temporarily in VaxPress/vaxpress/ directory
    tmp_outfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'LP_result.txt')
    if os.path.exists(tmp_outfile):
        os.remove(tmp_outfile)

    cmd = [linearpartition_bin, '-o', tmp_outfile]
    with sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE) as proc:
        proc.stdin.write(sequence.encode())
        proc.stdin.write(b'\n')
        proc.stdin.close()
        #for chunk in read_live_updates(proc.stdout):
        #    print(chunk, end='')
        #for chunk in read_live_updates(proc.stderr):
        #    print(chunk, end='')
        if (proc.wait() != 0):
            raise RuntimeError('linearpartition failed with exit code {}'.format(proc.returncode))

    return tmp_outfile

def parse_linearpartition_output(tmp_outfile, n):
    # initialize Pi
    # Pi: dictionary of base pairing probabilities : key = i : value = Pi (i = 1, ..., n)
    #n = len(sequence)
    Pi = {i: 0 for i in range(1, n + 1)} # _length_cds ??
    tmp_outfile = os.path.abspath(tmp_outfile)

    # read output
    with open(tmp_outfile, "r") as f:
        while True:
            line = f.readline().strip()  # line[0]: i, line[1]: j, line[2]: Pij
            if not line:
                break
            line = line.split(" ")

            # calculate Pi
            Pi[int(line[0])] += float(line[2])
            Pi[int(line[1])] += float(line[2])
    os.remove(tmp_outfile) #
    return Pi

def calculate_aup(Pi, sequence):
    # calculate loop metric = sum(1-Pi) / n (n=sequence length)
    n = len(sequence)
    total_unpaired_prob = sum(1 - Pi[i] for i in Pi)
    average_unpaired_probability = total_unpaired_prob/n

    # calculate unpaired_ucount
    ucount = sequence.count("U")
    unpaired_ucount = ucount * average_unpaired_probability

    return average_unpaired_probability, unpaired_ucount

def get_pairingprob(sequence):
    linearpartition_dir = os.path.abspath("/qbio/shjo/2_vaccine_design/base_pairing/LinearPartition") # later needs to be changed
    tmp_outfile = run_linearpartition(linearpartition_dir, sequence)
    Pi = parse_linearpartition_output(tmp_outfile, len(sequence))
    return calculate_aup(Pi, sequence)