#!/usr/bin/env python3

import argparse
import pyfaidx
import pyfastx
import re
import sys
import itertools
import operator
import collections
import Bio.Seq


parser = argparse.ArgumentParser(description='Generate consensus eccDNA sequences from rolling circle reads')

parser.add_argument('--fastq',      type=str, required=True, help='input reads in fastq format')
parser.add_argument('--paf',        type=str, required=True, help='input alignments in PAF format')
parser.add_argument('--info',       type=str, required=True, help='output file for sequences information')
parser.add_argument('--seq',        type=str, required=True, help='output file for consensus sequences')
parser.add_argument('--var',        type=str, required=True, help='output file for variants')
parser.add_argument('--reference',  type=str, required=True, help='reference genome sequences in fasta format')
parser.add_argument('--verbose',    action='store_true', help='print details of consensus construction')

parser.add_argument('--maxOffset',  type=int,   default=20,   help='maximum offset of start/end positions between two sub-reads to be considered as mapping to the same location')
parser.add_argument('--minMapQual', type=int,   default=30,   help='minimum mapping quality of sub-reads')
parser.add_argument('--minDP',      type=int,   default=4,    help='minimum depth to call variants')
parser.add_argument('--minAF',      type=float, default=0.75, help='minimum alternative allele frequency to call variants')

args = parser.parse_args()


def multimode(data):
    counts = collections.Counter(iter(data)).most_common()
    _, mode_items = next(itertools.groupby(counts, key=operator.itemgetter(1)), (0, []))
    return list(map(operator.itemgetter(0), mode_items))


class Interval(object):
    def __init__(self, chr, start, end, strand):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand


class Fragment(object):
    def __init__(self, start, end, chr, chrStart, chrEnd, strand, mapQual, cigar):
        # Both start and end are 0-based and inclusive
        self.start = start
        self.end = end
        self.interval = Interval(chr, chrStart, chrEnd, strand)
        self.mapQual = mapQual
        self.cigar = cigar


class FragmentGroup(object):
    def __init__(self, maxOffset=args.maxOffset):
        self.group = []
        self.intervals = []
        self.nfrag = 0
        self.passes = 0

        self.pointer = 0
        self.addFirst = False
        self.addLast = False

        self.maxOffset = maxOffset
    
    def _sameLocation(self, index, interval, withStart=True, withEnd=True):
        if len(self.group[index]) == 0:
            return True
        if interval.chr != self.group[index][0].interval.chr:
            return False
        if withStart:
            start = min(multimode([frag.interval.start for frag in self.group[index]]))
            if abs(start - interval.start) > self.maxOffset:
                return False
        if withEnd:
            end = max(multimode([frag.interval.end for frag in self.group[index]]))
            if abs(end - interval.end) > self.maxOffset:
                return False
        return True
    
    def _rewind(self, interval):
        return self._sameLocation(0, interval)
    
    def _noOverlap(self, interval):
        for frag in self.group:
            if frag[0].interval.chr == interval.chr:
                if interval.start <= frag[0].interval.end and interval.end >= frag[0].interval.start:
                    return False
        return True
    
    def add(self, fragment):
        if len(self.group) == 0:
            self.group.append([fragment])
            return True
        if self.nfrag == 0:
            if self._rewind(fragment.interval):
                self.group[0].append(fragment)
                self.nfrag = len(self.group)
                self.passes = 1
                if self.pointer == 0:
                    self.passes += 1
                self.pointer = 0
                return True
            else:
                if self._noOverlap(fragment.interval):
                    self.group.append([fragment])
                    self.pointer += 1
                    return True
        else:
            pointer = self.pointer + 1
            if pointer == self.nfrag:
                pointer = 0
            if self._sameLocation(pointer, fragment.interval):
                self.group[pointer].append(fragment)
                self.pointer = pointer
                if self.pointer == self.nfrag - 1:
                    self.passes += 1
                return True
        return False
    
    def addTerminals(self, firstFragment, lastFragment):
        if len(self.group) == 0:
            return
        self.addLast = False
        self.addFirst = False
        # Add last fragment
        if self.nfrag > 0:
            pointer = self.pointer + 1
            if pointer == self.nfrag:
                pointer = 0
        else:
            pointer = 0
        if lastFragment.interval.strand == '+':
            withStart = True
            withEnd = False
        else:
            withStart = False
            withEnd = True
        if self._sameLocation(pointer, lastFragment.interval, withStart, withEnd):
            self.addLast = True
            self.pointer = pointer
        # Add first fragment
        if firstFragment.interval.strand == '+':
            withStart = False
            withEnd = True
        else:
            withStart = True
            withEnd = False
        if self._sameLocation(-1, firstFragment.interval, withStart, withEnd):
            self.addFirst = True
        
        if self.addLast:
            self.group[pointer].append(lastFragment)
            if self.nfrag == 0:
                self.nfrag = len(self.group)
                self.passes = 1
        if self.addFirst:
            self.group[-1].append(firstFragment)
            if self.nfrag == 0:
                self.nfrag = len(self.group)
                self.passes = 1
        if self.addLast and self.addFirst:
            if pointer == self.nfrag - 1 and firstFragment.interval.strand == lastFragment.interval.strand:
                if firstFragment.interval.strand == '+':
                    if firstFragment.interval.start <= lastFragment.interval.end + 1:
                        self.passes += 1
                else:
                    if lastFragment.interval.start <= firstFragment.interval.end + 1:
                        self.passes += 1
    
    def reduce(self):
        if self.nfrag == 0:
            return None
        for frags in self.group:
            chr = frags[0].interval.chr
            start = min(multimode([frag.interval.start for frag in frags]))
            end = max(multimode([frag.interval.end for frag in frags]))
            strand = multimode([frag.interval.strand for frag in frags])[0]
            self.intervals.append(Interval(chr, start, end, strand))
        return self
    
    def print(self, read):
        sys.stdout.write(read.name + '\n')
        sys.stdout.write("#Fragment: " + str(self.nfrag) + "\tFull Pass: " + str(self.passes) + '\tRead Length: ' + str(read.length) + '\n\n')
        
        ntotal = 0
        for frag in self.group:
            ntotal += len(frag)
        if self.addFirst:
            ntotal -= 1
        
        if self.addFirst:
            frag = self.group[-1][-1]
            sys.stdout.write('\t' * self.nfrag)
            sys.stdout.write('{0: >5}'.format(frag.start+1) + ' - ' + frag.interval.chr + ':' + str(frag.interval.start+1) + '-' + str(frag.interval.end+1) + ' (' + frag.interval.strand + ') - ' + '{0: >5}'.format(frag.end+1))
            sys.stdout.write('\n')
        for i in range(1, ntotal+1):
            frag = self.group[(i-1) % self.nfrag][(i-1) // self.nfrag]
            sys.stdout.write('\t{0: >5}'.format(frag.start+1) + ' - ' + frag.interval.chr + ':' + str(frag.interval.start+1) + '-' + str(frag.interval.end+1) + ' (' + frag.interval.strand + ') - ' + '{0: >5}'.format(frag.end+1))
            if i % self.nfrag == 0:
                sys.stdout.write('\n')
        if ntotal % self.nfrag != 0:
            sys.stdout.write('\n')
        
        sys.stdout.write('\nLocation:\n')
        for interval in self.intervals:
            sys.stdout.write('\t' + interval.chr + ':' + str(interval.start+1) + '-' + str(interval.end+1) + ' (' + interval.strand + ')')
        sys.stdout.write('\n\n\n')


class Read(object):
    def __init__(self, name=None, length=0, minMapQual=args.minMapQual):
        self.name = name
        self.length = int(length)
        self.fragments = []
        self.minMapQual = minMapQual
    
    def isEmpty(self):
        return len(self.fragments) == 0
    
    def add(self, fragment):
        self.fragments.append(fragment)
    
    def ordered(self):
        self.fragments.sort(key=lambda x: x.start)
        return self
    
    def threadFix(self):
        self.ordered()
        for i in range(len(self.fragments)-1):
            dislen = self.fragments[i+1].start - self.fragments[i].end - 1
            if dislen > 0:
                if self.fragments[i].interval.strand == '+':
                    self.fragments[i].end += dislen
                    self.fragments[i].interval.end += dislen
                    self.fragments[i].cigar += '{0}M'.format(dislen)
                else:
                    if self.fragments[i+1].interval.strand == '-':
                        self.fragments[i+1].start -= dislen
                        self.fragments[i+1].interval.end += dislen
                        self.fragments[i+1].cigar += '{0}M'.format(dislen)
                    else:
                        self.fragments[i+1].start -= dislen
                        self.fragments[i+1].interval.start -= dislen
                        self.fragments[i+1].cigar = '{0}M'.format(dislen) + self.fragments[i+1].cigar
            elif dislen < 0:
                dislen = -dislen
                if self.fragments[i].interval.strand == '+':
                    self.fragments[i].end -= dislen
                    cigar_re = re.compile(r"(\d+)([MDI])")
                    ces = [[int(m.group(1)), m.group(2)] for m in cigar_re.finditer(self.fragments[i].cigar)]
                    p = len(ces) - 1
                    while dislen > 0 and p >= 0:
                        if ces[p][1] == 'M' or ces[p][1] == 'I':
                            s = min(ces[p][0], dislen)
                            dislen -= s
                            ces[p][0] -= s
                            if ces[p][1] == 'M':
                                self.fragments[i].interval.end -= s
                            if ces[p][0] == 0:
                                p -= 1
                        if p >= 0 and ces[p][1] == 'D':
                            self.fragments[i].interval.end -= ces[p][0]
                            ces[p][0] = 0
                            p -= 1
                    if p < 0:
                        return False
                    if ces[p][1] == 'I':
                        ces[p][1] = 'M'
                        self.fragments[i].interval.end += ces[p][0]
                    self.fragments[i].cigar = ''.join([str(step) + op for step, op in ces[0:p+1]])
                else:
                    if self.fragments[i+1].interval.strand == '-':
                        self.fragments[i+1].start += dislen
                        cigar_re = re.compile(r"(\d+)([MDI])")
                        ces = [[int(m.group(1)), m.group(2)] for m in cigar_re.finditer(self.fragments[i+1].cigar)]
                        p = len(ces) - 1
                        while dislen > 0 and p >= 0:
                            if ces[p][1] == 'M' or ces[p][1] == 'I':
                                s = min(ces[p][0], dislen)
                                dislen -= s
                                ces[p][0] -= s
                                if ces[p][1] == 'M':
                                    self.fragments[i+1].interval.end -= s
                                if ces[p][0] == 0:
                                    p -= 1
                            if p >= 0 and ces[p][1] == 'D':
                                self.fragments[i+1].interval.end -= ces[p][0]
                                ces[p][0] = 0
                                p -= 1
                        if p < 0:
                            return False
                        if ces[p][1] == 'I':
                            ces[p][1] = 'M'
                            self.fragments[i+1].interval.end += ces[p][0]
                        self.fragments[i+1].cigar = ''.join([str(step) + op for step, op in ces[0:p+1]])
                    else:
                        self.fragments[i+1].start += dislen
                        cigar_re = re.compile(r"(\d+)([MDI])")
                        ces = [[int(m.group(1)), m.group(2)] for m in cigar_re.finditer(self.fragments[i+1].cigar)]
                        p = 0
                        while dislen > 0 and p < len(ces):
                            if ces[p][1] == 'M' or ces[p][1] == 'I':
                                s = min(ces[p][0], dislen)
                                dislen -= s
                                ces[p][0] -= s
                                if ces[p][1] == 'M':
                                    self.fragments[i+1].interval.start += s
                                if ces[p][0] == 0:
                                    p += 1
                            if p < len(ces) and ces[p][1] == 'D':
                                self.fragments[i+1].interval.start += ces[p][0]
                                ces[p][0] = 0
                                p += 1
                        if p >= len(ces):
                            return False
                        if ces[p][1] == 'I':
                            ces[p][1] = 'M'
                            self.fragments[i+1].interval.start -= ces[p][0]
                        self.fragments[i+1].cigar = ''.join([str(step) + op for step, op in ces[p:]])
        return True
    
    def bootstrap(self):
        if len(self.fragments) < 3:
            return None
        res = self.threadFix()
        if res == False:
            return None
        fraggroup = FragmentGroup()
        for i in range(1, len(self.fragments)-1):
            if self.fragments[i].mapQual < self.minMapQual:
                break
            res = fraggroup.add(self.fragments[i])
            if res == False:
                break
        fraggroup.addTerminals(self.fragments[0], self.fragments[-1])
        fraggroup = fraggroup.reduce()
        return fraggroup


class PAF(object):
    def __init__(self, filepath):
        self.filepath = filepath
    
    def reads(self):
        curread = ''
        read = Read()
        with open(self.filepath, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if fields[0] != curread:
                    if not read.isEmpty():
                        yield read
                    curread = fields[0]
                    read = Read(name=fields[0], length=fields[1])
                assert fields[-1][0:2] == 'cg', 'The last column in PAF should be CIGAR string'
                frag = Fragment(start=int(fields[2]), end=int(fields[3])-1,
                                chr=fields[5], chrStart=int(fields[7]), chrEnd=int(fields[8])-1, strand=fields[4],
                                mapQual=int(fields[11]), cigar=fields[-1].split(':')[2])
                read.add(frag)
            if not read.isEmpty():
                yield read


class Base(object):
    def __init__(self):
        self.nt = {}
        self.last = ''
    
    def match(self, nt):
        if nt in self.nt:
            self.nt[nt] += 1
        else:
            self.nt[nt] = 1
        self.last = nt
    
    def insert(self, nt):
        self.nt[self.last] -= 1
        if self.last != '-':
            nt = self.last + nt
        if nt in self.nt:
            self.nt[nt] += 1
        else:
            self.nt[nt] = 1
    
    def call(self):
        alt = None
        cnt = 0
        depth = 0
        for nt in self.nt:
            depth += self.nt[nt]
            if self.nt[nt] > cnt:
                alt = nt
                cnt = self.nt[nt]
        return alt, cnt, depth


class Variant(object):
    def __init__(self, chr, pos, ref, alt, count, depth):
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.count = count
        self.depth = depth


class SequenceVariant(object):
    def __init__(self, reference, minAF=args.minAF, minDP=args.minDP):
        self.reference = reference
        self.sequence = [Base() for _ in range(len(reference))]
        
        self.minAF = minAF
        self.minDP = minDP
    
    def align(self, seq, readInterval, refInterval, cigar):
        if readInterval.chr != refInterval.chr or readInterval.start > refInterval.end or refInterval.start > readInterval.end:
            return False
        cigar_re = re.compile(r"(\d+)([MDI])")
        ces = [(int(m.group(1)), m.group(2)) for m in cigar_re.finditer(cigar)]
        read_cur = 0
        ref_cur = readInterval.start - refInterval.start
        for step, op in ces:
            if op == 'M':
                for _ in range(step):
                    if ref_cur >= 0 and ref_cur < len(self.sequence):
                        self.sequence[ref_cur].match(seq[read_cur])
                    read_cur += 1
                    ref_cur += 1
            elif op == 'D':
                for _ in range(step):
                    if ref_cur >= 0 and ref_cur < len(self.sequence):
                        self.sequence[ref_cur].match('-')
                    ref_cur += 1
            elif op == 'I':
                if ref_cur > 0 and ref_cur <= len(self.sequence):
                    self.sequence[ref_cur-1].insert(seq[read_cur:read_cur+step])
                read_cur += step
        return True
    
    def qualified(self, var, pos):
        if var.depth < self.minDP:
            return False
        if float(var.count) / var.depth < self.minAF:
            return False
        # Filter del and ins in homopolymer region
        homoins = False
        if len(var.alt) > 1:
            homoins = True
            for i in range(1, len(var.alt)):
                if var.alt[i] != var.ref:
                    homoins = False
                    break
        if var.alt == '-' or homoins:
            if pos+2 < len(self.reference) and self.reference[pos+1] == var.ref and self.reference[pos+2] == var.ref:
                return False
            if pos-1 >= 0 and pos+1 < len(self.reference) and self.reference[pos-1] == var.ref and self.reference[pos+1] == var.ref:
                return False
            if pos-2 >= 0 and self.reference[pos-1] == var.ref and self.reference[pos-2] == var.ref:
                return False
        return True
    
    def getSeqVar(self, refInterval):
        seq = [''] * len(self.reference)
        var = []
        for i in range(len(self.reference)):
            alt, cnt, depth = self.sequence[i].call()
            ref = self.reference[i]
            if alt is None or alt == ref:
                seq[i] = ref
            else:
                v = Variant(refInterval.chr, refInterval.start + i, ref, alt, cnt, depth)
                if self.qualified(v, i):
                    if alt != '-':
                        seq[i] = alt
                        var.append(v)
                    else:
                        if len(var) > 0 and var[-1].alt == '-':
                            var[-1].ref += ref
                            var[-1].count = min(var[-1].count, cnt)
                            var[-1].depth = min(var[-1].depth, depth)
                        else:
                            var.append(v)
                else:
                    seq[i] = ref
        seq = ''.join(seq)
        return seq, var


class ConsensusSequence(object):
    def __init__(self, fastq, genome, read):
        self.fastq = fastq
        self.genome = genome
        self.read = read
        self.fraggroup = read.bootstrap()
        
        self.sequence = ''
        self.variants = []
    
    def _getReference(self):
        seq = []
        if self.fraggroup is not None:
            for interval in self.fraggroup.intervals:
                s = self.genome[interval.chr][interval.start:interval.end+1]
                seq.append(s.seq)
        return seq
    
    def callConsensus(self):
        reference = self._getReference()
        if not reference:
            return self
        readseq = self.fastq[self.read.name].seq
        consensus = [SequenceVariant(seq) for seq in reference]
        
        for i in range(len(reference)):
            for frag in self.fraggroup.group[i]:
                seq = readseq[frag.start:frag.end+1]
                if frag.interval.strand == '-':
                    seq = str(Bio.Seq.Seq(seq).reverse_complement())
                readInterval = frag.interval
                refInterval = self.fraggroup.intervals[i]
                consensus[i].align(seq, readInterval, refInterval, frag.cigar)
        
        for i in range(len(reference)):
            seq, var = consensus[i].getSeqVar(refInterval=self.fraggroup.intervals[i])
            if self.fraggroup.intervals[i].strand == '-':
                seq = str(Bio.Seq.Seq(seq).reverse_complement())
            self.sequence += seq
            self.variants += var
        return self
    
    def getFasta(self):
        fa = ''
        if self.sequence:
            fa = '>' + self.read.name + '\n' + self.sequence + '\n'
        return fa
    
    def getVariants(self):
        for var in self.variants:
            yield '\t'.join([var.chr, str(var.pos+1), var.ref, var.alt, str(var.count), str(var.depth)]) + '\n'
    
    def getInfo(self):
        readname = self.read.name
        passes = 0
        nfrag = 0
        reflength = 0
        seqlength = len(self.sequence)
        frags = ''
        if self.fraggroup is not None:
            passes = self.fraggroup.passes
            nfrag = self.fraggroup.nfrag
            frags = []
            for interval in self.fraggroup.intervals:
                reflength += interval.end - interval.start + 1
                frags.append(interval.chr + ':' + str(interval.start+1) + '-' + str(interval.end+1) + '(' + interval.strand + ')')
            frags = '|'.join(frags)
        return '\t'.join([readname, str(passes), str(nfrag), str(reflength), str(seqlength), frags]) + '\n'
    
    def print(self):
        if self.fraggroup is not None:
            self.fraggroup.print(self.read)




if __name__ == "__main__":
    with open(args.info, 'w') as out, open(args.seq, 'w') as fa, open(args.var, 'w') as var:
        out.write('\t'.join(('readname', 'Nfullpass', 'Nfragment', 'refLength', 'seqLength', 'fragments')) + '\n')
        
        fastq = pyfastx.Fastq(args.fastq)
        genome = pyfaidx.Fasta(args.reference)
        
        for read in PAF(args.paf).reads():
            cs = ConsensusSequence(fastq, genome, read).callConsensus()
            seq = cs.getFasta()
            info = cs.getInfo()
            fa.write(seq)
            out.write(info)
            for v in cs.getVariants():
                var.write(v)
            if args.verbose:
                cs.print()



