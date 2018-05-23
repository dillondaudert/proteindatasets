# NOTE: This requires the DSSPData class shown in this gist: https://gist.github.com/jlhg/5181883
import os
import pandas as pd
from multiprocessing import Process, Queue
from DSSPData import DSSPData

datapath = "/home/dillon/data/dssp/raw/"

def _parse(filename):
    '''Parse a single .dssp file.
    
    This function splits a dssp record into its constituent
    chains and returns a record for each chain.
    
    Returns a list of records, each with:
        dssp_chain_id           - str
        seq                     - str
        ss                      - str
        acc                     - csv string
        tco                     - csv string
        kappa                   - csv string
        alpha                   - csv string
        phi                     - csv string
        psi                     - csv string
    '''

    dssp = DSSPData()
    dssp.parseDSSP(filename)
    dssp_id = os.path.splitext(os.path.basename(filename))[0]
    
    # get lists of features for entire entry
    seq = dssp.aa
    structs = []
    for i in range(len(dssp.aa)):
        s = dssp.struct[i][2]
        if s == ' ':
            if dssp.struct[i][0] == "*":
                assert seq[i] == "!"
                seq[i] = "*"
            s = 'U'
        structs.append(s)
    ss = structs
    acc = dssp.acc
    tco = dssp.tco
    kappa = dssp.kappa
    alpha = dssp.alpha
    phi = dssp.phi
    psi = dssp.psi
    
    # check lengths
    try:
        assert len(seq) == len(structs) == len(acc) == len(tco) 
        assert len(seq) == len(kappa) == len(alpha) == len(phi)
        assert len(seq) == len(psi)
    except AssertionError as e:
        print("In dssp record %s, field lengths did not match.\n" % dssp_id)
        raise
        
    chain_type = list(dssp.getChainType())
    chain_seps = []
    for i, aa in enumerate(seq):
        # find the start and end locations
        if aa == "*":
            chain_seps.append(i)
            
    chain_names = set(chain_type)
    if "" in chain_names:
        chain_names.remove("")
        
    try:
        assert "".join(list(chain_names)).isalnum()
        assert len(chain_names) == len(chain_seps)+1
    except AssertionError:
        print(dssp_id, "chain_names: ", chain_names, "chain_seps: ", chain_seps)
        raise
    
    chain_records = []
    unique_chain_ids = set()
    
    prev_chain_end = -1    
    for i, chain_end in enumerate(chain_seps):
        chain_id = chain_type[prev_chain_end+1]
        chain_dssp_id = dssp_id + chain_id
        unique_chain_ids.add(chain_dssp_id)
        chain_seq = "".join(seq[prev_chain_end+1:chain_end])
        chain_ss = "".join(ss[prev_chain_end+1:chain_end])
        chain_acc = ",".join(acc[prev_chain_end+1:chain_end])
        chain_tco = ",".join(tco[prev_chain_end+1:chain_end])
        chain_kappa = ",".join(kappa[prev_chain_end+1:chain_end])
        chain_alpha = ",".join(alpha[prev_chain_end+1:chain_end])
        chain_phi = ",".join(phi[prev_chain_end+1:chain_end])
        chain_psi = ",".join(psi[prev_chain_end+1:chain_end])
        chain_records.append((chain_dssp_id, chain_seq, chain_ss, 
                              chain_acc, chain_tco, chain_kappa, 
                              chain_alpha, chain_phi, chain_psi))
        prev_chain_end = chain_end
        
    # then the last chain
    chain_dssp_id = dssp_id + chain_type[prev_chain_end+1]
    unique_chain_ids.add(chain_dssp_id)
    chain_seq = "".join(seq[prev_chain_end+1:])
    chain_ss = "".join(ss[prev_chain_end+1:])
    chain_acc = ",".join(acc[prev_chain_end+1:])
    chain_tco = ",".join(tco[prev_chain_end+1:])
    chain_kappa = ",".join(kappa[prev_chain_end+1:])
    chain_alpha = ",".join(alpha[prev_chain_end+1:])
    chain_phi = ",".join(phi[prev_chain_end+1:])
    chain_psi = ",".join(psi[prev_chain_end+1:])
    chain_records.append((chain_dssp_id, chain_seq, chain_ss, 
                          chain_acc, chain_tco, chain_kappa, 
                          chain_alpha, chain_phi, chain_psi))
    
    try:
        assert len(chain_names) == len(unique_chain_ids)
    except AssertionError:
        print("chain names: ", chain_names, "\nunique ids: ", unique_chain_ids)
        raise
        
    
    return chain_records
        
        
        
        


def parser(worker_queue, done_queue, id):
    '''Parser task that parses DSSP files from the worker_queue and places
    the results in the done_queue.
    '''

    local_count = 0

    while True:
        filename = worker_queue.get()
        # check if done
        if filename is None:
            return

        local_count += 1
        if local_count % 1000 == 0:
            print("Process #%d has handled %d records." % (id, local_count))

        try:
            records = _parse(filename)
            done_queue.put(records)
        except:
            print("Worker %d encountered exception, skipping record...")


if __name__ == '__main__':

    num_workers = 4
    files = os.listdir(datapath)
    num_records = len(files)

    print("Processing %d individual records." % (num_records))

    worker_queue = Queue()
    done_queue = Queue()

    print("Controller's pid is %d" % (os.getpid()))
    print("Spawning %d workers..." % (num_workers))
    workers = []
    for i in range(num_workers):
        p = Process(target=parser, args=(worker_queue, done_queue, i))
        workers.append(p)
        p.start()

    # fill the queue with files
    for filename in files:
        worker_queue.put(datapath+filename)

    # Add sentinels to signal end
    for _ in range(num_workers):
        worker_queue.put(None)

    print("Parsing DSSP files...")
    results = []
    for i in range(num_records):
        records = done_queue.get()
        for rec in records:
            results.append(rec)

    assert len(results) != 0

    print("Saving %d records to dssp_records.csv" % len(results))
    filename = "dssp_records.csv"
    df = pd.DataFrame.from_records(results, columns=["dssp_id", "seq", "ss", "acc", "tco", "kappa", "alpha", "phi", "psi"])
    df.to_csv(filename, index=False)

    print("Joining %d workers..." % (num_workers))
    for p in workers:
        p.join()
