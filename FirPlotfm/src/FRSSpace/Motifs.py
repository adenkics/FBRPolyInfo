import csv
import pickle, os, copy, re, random

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem import rdmolops
from rdkit.Chem import PandasTools
from FirPlotfm.src.utils.data_utils import VCGData


script_dir = os.path.dirname(os.path.abspath(__file__))
relative_path = os.path.join('..', '..', 'data', 'MotifDefsData.pkl')
pkl_file_path = os.path.abspath(os.path.join(script_dir, relative_path))

with open(pkl_file_path, 'rb') as file:
  MotifDefs = pickle.load(file)

environs = MotifDefs['environs']
reactionDefs = MotifDefs['reactionDefs']

smartsGps = copy.deepcopy(reactionDefs)

for gp in smartsGps:
    for j, defn in enumerate(gp):
        g1, g2, bnd = defn
        r1 = environs['L' + g1]
        r2 = environs['L' + g2]
        g1 = re.sub('[a-z,A-Z]', '', g1)
        g2 = re.sub('[a-z,A-Z]', '', g2)
        sma = '[$(%s):1]%s;!@[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]' % (r1, bnd, r2, g1, g2)
        gp[j] = sma

for gp in smartsGps:
    for defn in gp:
        try:
            t = Reactions.ReactionFromSmarts(defn)
            t.Initialize()
        except Exception:
            print(defn)
            raise

environMatchers = {}
for env, sma in environs.items():
    environMatchers[env] = Chem.MolFromSmarts(sma)

bondMatchers = []
for i, compats in enumerate(reactionDefs):
    tmp = []
    for i1, i2, bType in compats:
        e1 = environs['L%s' % i1]
        e2 = environs['L%s' % i2]
        patt = '[$(%s)]%s;!@[$(%s)]' % (e1, bType, e2)
        patt = Chem.MolFromSmarts(patt)
        tmp.append((i1, i2, bType, patt))
    bondMatchers.append(tmp)

reactions = tuple([[Reactions.ReactionFromSmarts(y) for y in x] for x in smartsGps])
reverseReactions = []
for i, rxnSet in enumerate(smartsGps):
    for j, sma in enumerate(rxnSet):
        rs, ps = sma.split('>>')
        sma = '%s>>%s' % (ps, rs)
        rxn = Reactions.ReactionFromSmarts(sma)
        labels = re.findall(r'\[([0-9]+?)\*\]', ps)
        rxn._matchers = [Chem.MolFromSmiles('[%s*]' % x) for x in labels]
        reverseReactions.append(rxn)

def MotifDecomp(mol, allNodes=None, minFragmentSize=1, onlyUseReactions=None, silent=True,
                keepNonLeafNodes=False, singlePass=False, returnMols=False):
  """ returns the Motifs decomposition for a pyrolysis prodcuts
  """
  global reactions
  mSmi = Chem.MolToSmiles(mol, 1)

  if allNodes is None:
    allNodes = set()

  if mSmi in allNodes:
    return set()

  activePool = {mSmi: mol}
  allNodes.add(mSmi)
  foundMols = {mSmi: mol}
  for gpIdx, reactionGp in enumerate(reactions):
    newPool = {}
    while activePool:
      matched = False
      nSmi = next(iter(activePool))
      mol = activePool.pop(nSmi)
      for rxnIdx, reaction in enumerate(reactionGp):
        if onlyUseReactions and (gpIdx, rxnIdx) not in onlyUseReactions:
          continue
        if not silent:
          print('--------')
          print(smartsGps[gpIdx][rxnIdx])
        ps = reaction.RunReactants((mol,))
        if ps:
          if not silent:
            print(nSmi, '->', len(ps), 'products')
          for prodSeq in ps:
            seqOk = True
            # we want to disqualify small fragments, so sort the product sequence by size
            tSeq = [(prod.GetNumAtoms(onlyExplicit=True), idx)
                    for idx, prod in enumerate(prodSeq)]
            tSeq.sort()
            for nats, idx in tSeq:
              prod = prodSeq[idx]
              try:
                Chem.SanitizeMol(prod)
              except Exception:
                continue
              pSmi = Chem.MolToSmiles(prod, 1)
              if minFragmentSize > 0:
                nDummies = pSmi.count('*')
                if nats - nDummies < minFragmentSize:
                  seqOk = False
                  break
              prod.pSmi = pSmi
            ts = [(x, prodSeq[y]) for x, y in tSeq]
            prodSeq = ts
            if seqOk:
              matched = True
              for nats, prod in prodSeq:
                pSmi = prod.pSmi
                # print('\t',nats,pSmi)
                if pSmi not in allNodes:
                  if not singlePass:
                    activePool[pSmi] = prod
                  allNodes.add(pSmi)
                  foundMols[pSmi] = prod
      if singlePass or keepNonLeafNodes or not matched:
        newPool[nSmi] = mol
    activePool = newPool
  if not (singlePass or keepNonLeafNodes):
    if not returnMols:
      res = set(activePool.keys())
    else:
      res = activePool.values()
  else:
    if not returnMols:
      res = allNodes
    else:
      res = foundMols.values()
  return res
