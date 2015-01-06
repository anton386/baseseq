import sys


class Chain(object):

    def __init__(self, chain):
        self.chain = chain

    def output_chain(self, ref, consensus, structure):

        #TODO raise exception if reference and consensus are not the same
        
        f_out = open(self.chain, "w")
        
        identifier = "chain"
        score = 5000
        tname = "reference"
        tsize = len(ref.sequence)
        tstrand = "+"
        tstart = 0
        tend = len(ref.sequence)
        qname = "consensus"
        qsize = len(consensus)
        qstrand = "+"
        qstart = 0
        qend = len(consensus)
        chain_id = 1

        f_out.write("\t".join(map(str,
                                  [identifier, score,
                                   tname, tsize,
                                   tstrand, tstart, tend,
                                   qname, qsize,
                                   qstrand, qstart, qend,
                                   chain_id])) + "\n")

        size = 0
        dt = 0  # diff btw end of block to next block (ref)
        dq = 0  # diff btw end of block to next block (query)
        for pos, (r, q) in enumerate(zip(ref.sequence, structure["sequence"])):
            position = pos + 1
            
            if q == "D":
                dq += 1
            elif q == "I":
                # dt += indels[position-1][1]
                dt += (len(structure["ins"][position])-1)
            else:
                size += 1

                if dq > 0:
                    # write to line
                    f_out.write("%s\t%s\t%s\n" % (size, dt, dq))
                    dq = 0
                    dt = 0
                    size = 0

                if dt > 0:
                    # write to line
                    f_out.write("%s\t%s\t%s\n" % (size, dt, dq))
                    dq = 0
                    dt = 0
                    size = 0

        if dt == dq == 0:
            f_out.write("%s\n" % (size))

        f_out.close()
