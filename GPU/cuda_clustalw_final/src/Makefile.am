## Makefile.am
## ------------------------------------------------------
## Este Makefile.am supõe que você queira usar OpenMP via NVCC
## (passando -fopenmp ao compilador host g++ usando -Xcompiler).
## Se não quiser, basta remover as linhas marcadas com "OPENMP".

AUTOMAKE_OPTIONS = foreign subdir-objects

SUBDIRS =

# Recomendado substituir INCLUDES por AM_CPPFLAGS
AM_CPPFLAGS = -I@top_srcdir@/

# Caso queira remover OpenMP, remova/ajuste estas duas linhas:
NVCCFLAGS = -Xcompiler -fopenmp -g -G
OPENMP_LIBS = -lgomp

# Aqui vão as bibliotecas extras de link
# (CUDA + OpenMP)
clustalw2_LDADD = -L/usr/local/cuda/lib64 -lcuda -lcudart $(OPENMP_LIBS)

# Ajuste final de link. Observação:
# - Não coloque $(CXXFLAGS) aqui para não injetar -fopenmp "puro" no nvcc
clustalw2_LINK = $(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o $@

bin_PROGRAMS = clustalw2

# Compilando .cu e .o via nvcc
.cu.o:
	$(NVCC) -DHAVE_CONFIG_H -I. -I.. -o $@ -c $< $(NVCCFLAGS)

.o:
	$(NVCC) -DHAVE_CONFIG_H -I. -I.. -o $@ -c $< $(NVCCFLAGS)

# Fontes
clustalw2_SOURCES = \
  ./pairwise/FullPairwiseAlign.cpp \
  ./pairwise/FastPairwiseAlign.cpp \
  ./cuda/Stack.cu \
  ./cuda/cudaFullPairwiseAlign.cu \
  ./fileInput/MSFFileParser.cpp \
  ./fileInput/FileReader.cpp \
  ./fileInput/PIRFileParser.cpp \
  ./fileInput/RSFFileParser.cpp \
  ./fileInput/GDEFileParser.cpp \
  ./fileInput/InFileStream.cpp \
  ./fileInput/ClustalFileParser.cpp \
  ./fileInput/PearsonFileParser.cpp \
  ./fileInput/FileParser.cpp \
  ./fileInput/EMBLFileParser.cpp \
  ./tree/UPGMA/RootedClusterTree.cpp \
  ./tree/UPGMA/UPGMAAlgorithm.cpp \
  ./tree/UPGMA/Node.cpp \
  ./tree/UPGMA/RootedGuideTree.cpp \
  ./tree/UPGMA/RootedTreeOutput.cpp \
  ./tree/Tree.cpp \
  ./tree/ClusterTree.cpp \
  ./tree/TreeInterface.cpp \
  ./tree/UnRootedClusterTree.cpp \
  ./tree/ClusterTreeOutput.cpp \
  ./tree/RandomGenerator.cpp \
  ./tree/NJTree.cpp \
  ./tree/AlignmentSteps.cpp \
  ./interface/InteractiveMenu.cpp \
  ./interface/CommandLineParser.cpp \
  ./substitutionMatrix/SubMatrix.cpp \
  ./multipleAlign/Iteration.cpp \
  ./multipleAlign/MSA.cpp \
  ./multipleAlign/MyersMillerProfileAlign.cpp \
  ./multipleAlign/ProfileStandard.cpp \
  ./multipleAlign/ProfileWithSub.cpp \
  ./multipleAlign/ProfileBase.cpp \
  ./multipleAlign/LowScoreSegProfile.cpp \
  ./cuda/cudaMultiSeqAlign.cpp \
  ./general/OutputFile.cpp \
  ./general/UserParameters.cpp \
  ./general/Utility.cpp \
  ./general/InvalidCombination.cpp \
  ./general/DebugLog.cpp \
  ./general/ClustalWResources.cpp \
  ./general/VectorOutOfRange.cpp \
  ./general/SymMatrix.cpp \
  ./general/Stats.cpp \
  ./Help.cpp \
  ./alignment/Alignment.cpp \
  ./alignment/AlignmentOutput.cpp \
  ./alignment/ObjectiveScore.cpp \
  ./alignment/Sequence.cpp \
  ./Clustal.cpp \
  ./main.cpp \
  \
  ## E também os headers (não geram objetos, mas é bom listá-los).
  ./pairwise/PairwiseAlignBase.h \
  ./pairwise/FastPairwiseAlign.h \
  ./pairwise/FullPairwiseAlign.h \
  ./cuda/DyArray1D.h \
  ./cuda/DyArray2D.h \
  ./cuda/pairwiseAlignInfo.h \
  ./cuda/cudaFullPairwiseAlign.h \
  ./fileInput/FileParser.h \
  ./fileInput/EMBLFileParser.h \
  ./fileInput/MSFFileParser.h \
  ./fileInput/FileReader.h \
  ./fileInput/PIRFileParser.h \
  ./fileInput/RSFFileParser.h \
  ./fileInput/GDEFileParser.h \
  ./fileInput/InFileStream.h \
  ./fileInput/ClustalFileParser.h \
  ./fileInput/PearsonFileParser.h \
  ./tree/UPGMA/RootedClusterTree.h \
  ./tree/UPGMA/RootedGuideTree.h \
  ./tree/UPGMA/RootedTreeOutput.h \
  ./tree/UPGMA/UPGMAAlgorithm.h \
  ./tree/UPGMA/upgmadata.h \
  ./tree/UPGMA/Node.h \
  ./tree/Tree.h \
  ./tree/ClusterTree.h \
  ./tree/ClusterTreeAlgorithm.h \
  ./tree/TreeInterface.h \
  ./tree/dayhoff.h \
  ./tree/UnRootedClusterTree.h \
  ./tree/ClusterTreeOutput.h \
  ./tree/RandomGenerator.h \
  ./tree/NJTree.h \
  ./interface/InteractiveMenu.h \
  ./interface/CommandLineParser.h \
  ./substitutionMatrix/SubMatrix.h \
  ./substitutionMatrix/globalmatrix.h \
  ./substitutionMatrix/matrices.h \
  ./multipleAlign/ProfileBase.h \
  ./multipleAlign/LowScoreSegProfile.h \
  ./multipleAlign/Iteration.h \
  ./multipleAlign/MSA.h \
  ./multipleAlign/MyersMillerProfileAlign.h \
  ./multipleAlign/ProfileStandard.h \
  ./multipleAlign/ProfileWithSub.h \
  ./multipleAlign/ProfileAlignAlgorithm.h \
  ./cuda/cudaMultiSeqAlign.h \
  ./general/ClustalWResources.h \
  ./general/VectorOutOfRange.h \
  ./general/SymMatrix.h \
  ./general/clustalw.h \
  ./general/OutputFile.h \
  ./general/param.h \
  ./general/utils.h \
  ./general/UserParameters.h \
  ./general/SquareMat.h \
  ./general/Utility.h \
  ./general/SequenceNotFoundException.h \
  ./general/Array2D.h \
  ./general/RandomAccessLList.h \
  ./general/general.h \
  ./general/debuglogObject.h \
  ./general/statsObject.h \
  ./general/userparams.h \
  ./general/DebugLog.h \
  ./general/Stats.h \
  ./Help.h \
  ./alignment/Alignment.h \
  ./alignment/AlignmentOutput.h \
  ./alignment/ObjectiveScore.h \
  ./alignment/Sequence.h \
  ./Clustal.h \
  ./clustalw_version.h
