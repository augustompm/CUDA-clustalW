/**
 * SubMatrix.cpp (versão unificada para CPU ClustalW 2.0.11)
 * Adaptado da versão “GPU SubMatrix.cpp” + correções para CPU
 * 
 * Copyright (c) 2007-2009 Des Higgins, Julie Thompson, Toby Gibson
 * Licença GPL v2 (como ClustalW 2.x)
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <exception>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctype.h>
#include "SubMatrix.h"
#include "matrices.h"
#include "../general/InvalidCombination.cpp"
#include "../general/debuglogObject.h"

namespace clustalw
{
using namespace std;

/**
 * Construtor padrão
 */
SubMatrix::SubMatrix()
: sizenAAMatrix(276),
  sizeDNAMatrix(153),
  matrixAvgScore(0),
  QTDNAHistMatNum(DNAIUB),
  QTAAHistMatNum(AAHISTGONNETPAM250),
  QTsegmentDNAMatNum(DNAIUB),
  QTsegmentAAMatNum(QTAASEGGONNETPAM250)
{
    userSeries = false;
    setUpCrossReferences();

#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        logObject->logMsg("Creating the SubMatrix object\n");
    }
#endif

    try
    {
        // Alocar as matrizes “padronizadas”:
        blosum30mtVec  = new Matrix(blosum30mt,  blosum30mt  + sizenAAMatrix);
        blosum40mtVec  = new Matrix(blosum40mt,  blosum40mt  + sizenAAMatrix);
        blosum45mtVec  = new Matrix(blosum45mt,  blosum45mt  + sizenAAMatrix);
        blosum62mt2Vec = new Matrix(blosum62mt2, blosum62mt2 + sizenAAMatrix);
        blosum80mtVec  = new Matrix(blosum80mt,  blosum80mt  + sizenAAMatrix);

        pam20mtVec  = new Matrix(pam20mt,  pam20mt  + sizenAAMatrix);
        pam60mtVec  = new Matrix(pam60mt,  pam60mt  + sizenAAMatrix);
        pam120mtVec = new Matrix(pam120mt, pam120mt + sizenAAMatrix);
        pam350mtVec = new Matrix(pam350mt, pam350mt + sizenAAMatrix);

        idmatVec    = new Matrix(idmat, idmat + sizenAAMatrix);

        gon40mtVec  = new Matrix(gon40mt,  gon40mt  + sizenAAMatrix);
        gon80mtVec  = new Matrix(gon80mt,  gon80mt  + sizenAAMatrix);
        gon120mtVec = new Matrix(gon120mt, gon120mt + sizenAAMatrix);
        gon160mtVec = new Matrix(gon160mt, gon160mt + sizenAAMatrix);
        gon250mtVec = new Matrix(gon250mt, gon250mt + sizenAAMatrix);
        gon350mtVec = new Matrix(gon350mt, gon350mt + sizenAAMatrix);

        clustalvdnamtVec = new Matrix(clustalvdnamt, clustalvdnamt + sizeDNAMatrix);
        swgapdnamtVec    = new Matrix(swgapdnamt,    swgapdnamt    + sizeDNAMatrix);

        // user-defined e matrizes QT:
        userMat.resize(NUMRES*NUMRES);
        pwUserMat.resize(NUMRES*NUMRES);
        userDNAMat.resize(NUMRES*NUMRES);
        pwUserDNAMat.resize(NUMRES*NUMRES);
        QTscoreUserMatrix.resize(NUMRES*NUMRES);
        QTscoreUserDNAMatrix.resize(NUMRES*NUMRES);
        QTsegmentDNAMatrix.resize(NUMRES*NUMRES);
        QTsegmentAAMatrix.resize(NUMRES*NUMRES);

        // series
        userMatSeries.resize(MAXMAT);
        {
            auto fm = userMatSeries.begin();
            auto lm = userMatSeries.end();
            while(fm != lm)
            {
                fm->resize(NUMRES*NUMRES);
                ++fm;
            }
        }

        // cross ref series
        AAXrefseries.resize(MAXMAT);
        {
            auto fx = AAXrefseries.begin();
            auto lx = AAXrefseries.end();
            while(fx != lx)
            {
                fx->resize(NUMRES+1);
                ++fx;
            }
        }

        // Ajustar nomes default
        matrixNum     = 3;  // gonnet
        matrixName    = new string("gonnet");

        DNAMatrixNum  = 1;  // iub
        DNAMatrixName = new string("iub");

        pwMatrixNum   = 3;  // gonnet
        pwMatrixName  = new string("gonnet");

        pwDNAMatrixNum  = 1; // iub
        pwDNAMatrixName = new string("iub");
    }
    catch(const exception &ex)
    {
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
    }
}

SubMatrix::~SubMatrix()
{
    // Deletar os ponteiros para as matrizes “new”
    delete blosum30mtVec;
    delete blosum40mtVec;
    delete blosum45mtVec;
    delete blosum62mt2Vec;
    delete blosum80mtVec;
    delete pam20mtVec;
    delete pam60mtVec;
    delete pam120mtVec;
    delete pam350mtVec;
    delete idmatVec;
    delete gon40mtVec;
    delete gon80mtVec;
    delete gon120mtVec;
    delete gon160mtVec;
    delete gon250mtVec;
    delete gon350mtVec;
    delete clustalvdnamtVec;
    delete swgapdnamtVec;

    // deletar as strings
    delete matrixName;
    delete DNAMatrixName;
    delete pwMatrixName;
    delete pwDNAMatrixName;
}

/**
 * setValuesToDefault(): Restaura defaults
 */
void SubMatrix::setValuesToDefault()
{
    matrixAvgScore     = 0;
    QTDNAHistMatNum    = DNAIUB;
    QTAAHistMatNum     = AAHISTGONNETPAM250;
    QTsegmentDNAMatNum = DNAIUB;
    QTsegmentAAMatNum  = QTAASEGGONNETPAM250;
    userSeries         = false;

    // defaults
    matrixNum     = 3; // gonnet
    DNAMatrixNum  = 1; // iub
    pwMatrixNum   = 3; // gonnet
    pwDNAMatrixNum= 1; // iub
}

/**
 * setUpCrossReferences() – inicializa defaultAAXref e defaultDNAXref
 */
void SubMatrix::setUpCrossReferences()
{
    char c1, c2;
    short i, j, maxRes;

    // Tabela de aminoácidos de "matrices.h"
    string aminoAcidOrder   = "ABCDEFGHIKLMNPQRSTVWXYZ";
    // Tabela de nucleotídeos (IUB)
    string nucleicAcidOrder = "ABCDGHKMNRSTUVWXY";

    defaultAAXref.resize(NUMRES+1);
    defaultDNAXref.resize(NUMRES+1);

    // user-defined xrefs
    DNAXref.resize(NUMRES+1);
    AAXref.resize(NUMRES+1);
    pwAAXref.resize(NUMRES+1);
    pwDNAXref.resize(NUMRES+1);
    QTscoreXref.resize(NUMRES+1);
    QTscoreDNAXref.resize(NUMRES+1);
    QTsegmentDNAXref.resize(NUMRES+1);
    QTsegmentAAXref.resize(NUMRES+1);

    for(i=0; i<NUMRES; i++)
    {
        defaultAAXref[i]=-1;
        defaultDNAXref[i]=-1;
    }

    // crossref p/ amino:
    maxRes=0;
    for(i=0; (c1=aminoAcidOrder[i]); i++)
    {
        for(j=0; (c2=userParameters->getAminoAcidCode(j)); j++)
        {
            if(c1==c2)
            {
                defaultAAXref[i] = j;
                maxRes++;
                break;
            }
        }
        if((defaultAAXref[i]==-1) && (c1!='*'))
        {
            utilityObject->error("residue %c in matrices.h is not recognised", c1);
        }
    }

    // crossref p/ DNA:
    maxRes=0;
    for(i=0; (c1=nucleicAcidOrder[i]); i++)
    {
        for(j=0; (c2=userParameters->getAminoAcidCode(j)); j++)
        {
            if(c1==c2)
            {
                defaultDNAXref[i] = j;
                maxRes++;
                break;
            }
        }
        if((defaultDNAXref[i]==-1) && (c1!='*'))
        {
            utilityObject->error("nucleic acid %c in matrices.h is not recognised", c1);
        }
    }
}

/**
 * getPairwiseMatrix() - obtém a matriz a ser usada em pairwise alignment
 */
int SubMatrix::getPairwiseMatrix(int matrix[NUMRES][NUMRES],
                                 PairScaleValues &scale,
                                 int &matAvg)
{
    int _maxRes;
    Matrix *_matPtr;
    Xref   *_matXref;

#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        logObject->logMsg("In the function getPairwiseMatrix:\n");
    }
#endif

#ifdef OS_MAC
    scale.intScale = 10;
#else
    scale.intScale = 100;
#endif
    scale.gapOpenScale = 1.0f;
    scale.gapExtendScale=1.0f;

    if(userParameters->getDNAFlag())
    {
#if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            string msg = "    (DNA AND Pairwise) " + *pwDNAMatrixName + "\n";
            logObject->logMsg(msg);
        }
#endif
        // "iub" ou "clustalw" ou user
        if(pwDNAMatrixName->compare("iub")==0)
        {
            _matPtr = swgapdnamtVec;
            _matXref= &defaultDNAXref;
        }
        else if(pwDNAMatrixName->compare("clustalw")==0)
        {
            _matPtr = clustalvdnamtVec;
            _matXref= &defaultDNAXref;
            scale.gapOpenScale   = 0.6667f;
            scale.gapExtendScale = 0.751f;
        }
        else
        {
            // user-defined
            _matPtr = &pwUserDNAMat;
            _matXref= &pwDNAXref;
        }

        _maxRes = getMatrix(_matPtr, _matXref, matrix, true, scale.intScale);
        if(_maxRes==0) return -1;

        float _transitionWeight = userParameters->getTransitionWeight();
        // A=0, G=4, C=2, T=11 e 12... ajusta transitions
        matrix[0][4]   = (int)(_transitionWeight*matrix[0][0]);
        matrix[4][0]   = (int)(_transitionWeight*matrix[0][0]);
        matrix[2][11]  = (int)(_transitionWeight*matrix[0][0]);
        matrix[11][2]  = (int)(_transitionWeight*matrix[0][0]);
        matrix[2][12]  = (int)(_transitionWeight*matrix[0][0]);
        matrix[12][2]  = (int)(_transitionWeight*matrix[0][0]);
    }
    else
    {
#if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            string msg = "    (Protein AND Pairwise) " + *pwMatrixName + "\n";
            logObject->logMsg(msg);
        }
#endif

        // Se “blosum”, “pam”, “gonnet”, “id” ou user
        if(pwMatrixName->compare("blosum")==0)
        {
            // BLOSUM62
            _matPtr = blosum62mt2Vec;
            _matXref= &defaultAAXref;
        }
        else if(pwMatrixName->compare("blosum62")==0)
        {
            _matPtr = blosum62mt2Vec;
            _matXref= &defaultAAXref;
        }
        else if(pwMatrixName->compare("blosum45")==0)
        {
            _matPtr = blosum45mtVec;
            _matXref= &defaultAAXref;
        }
        else if(pwMatrixName->compare("pam")==0)
        {
            // pam350
            _matPtr = pam350mtVec;
            _matXref= &defaultAAXref;
        }
        else if(pwMatrixName->compare("gonnet")==0)
        {
            _matPtr = gon250mtVec;
            _matXref= &defaultAAXref;
            scale.intScale/=10;
        }
        else if(pwMatrixName->compare("id")==0)
        {
            _matPtr = idmatVec;
            _matXref= &defaultAAXref;
        }
        else
        {
            // user-defined
            _matPtr = &pwUserMat;
            _matXref= &pwAAXref;
        }

        _maxRes = getMatrix(_matPtr, _matXref, matrix, true, scale.intScale);
        if(_maxRes==0) return -1;
    }

#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        ostringstream outs;
        outs <<"    Called getMatrix.\n"
             <<"    intScale="<< scale.intScale
             <<", gapOpenScale="<< scale.gapOpenScale
             <<", gapExtendScale="<< scale.gapExtendScale<<"\n\n";
        logObject->logMsg(outs.str());
    }
#endif

    matAvg = matrixAvgScore;
    return _maxRes;
}

/**
 * getArgs(...) — **apenas uma** versão, evitando duplicatas
 */
int SubMatrix::getArgs(char *inline1, char *args[], int max)
{
    char *inptr = inline1;
    int i;
    for(i = 0; i < max; i++)
    {
        if((args[i] = strtok(inptr, " \t\n")) == NULL)
            break;
        inptr = NULL;
    }
    return i;
}

/**
 * getProfileAlignMatrix() – para Multiple Align
 */
int SubMatrix::getProfileAlignMatrix(int matrix[NUMRES][NUMRES],
                                     double pcid, int minLen,
                                     PrfScaleValues &scaleParam,
                                     int &matAvg)
{
    bool _negMatrix = userParameters->getUseNegMatrix();
    int maxRes=0;

    scaleParam.intScale=100;

#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        logObject->logMsg("In getProfileAlignMatrix:\n");
    }
#endif

    if(userParameters->getDNAFlag())
    {
#if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            ostringstream outs;
            outs<<"    (DNA AND Multiple align) "<< DNAMatrixName->c_str() <<"\n";
            logObject->logMsg(outs.str());
        }
#endif
        scaleParam.scale=1.0f;
        if(DNAMatrixName->compare("iub")==0)
        {
            _matPtr = swgapdnamtVec;
            _matXref= &defaultDNAXref;
        }
        else if(DNAMatrixName->compare("clustalw")==0)
        {
            _matPtr = clustalvdnamtVec;
            _matXref= &defaultDNAXref;
            scaleParam.scale=0.66f;
        }
        else
        {
            _matPtr = &userDNAMat;
            _matXref= &DNAXref;
        }
        maxRes = getMatrix(_matPtr, _matXref, matrix, _negMatrix, (int)scaleParam.intScale);
        if(maxRes==0) return -1;

        float _transitionWeight=userParameters->getTransitionWeight();
        matrix[(*_matXref)[0]][(*_matXref)[4]]  = (int)(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[4]][(*_matXref)[0]]  = (int)(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[2]][(*_matXref)[11]] = (int)(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[11]][(*_matXref)[2]] = (int)(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[2]][(*_matXref)[12]] = (int)(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[12]][(*_matXref)[2]] = (int)(_transitionWeight*matrix[0][0]);
    }
    else
    {
#if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            ostringstream outs;
            outs<<"    (Protein AND Multiple align) "<< matrixName->c_str()<<"\n";
            logObject->logMsg(outs.str());
        }
#endif
        scaleParam.scale=0.75f;
        if(matrixName->compare("blosum")==0)
        {
            if(_negMatrix || !userParameters->getDistanceTree())
            {
                _matPtr=blosum40mtVec;
            }
            else if(pcid>80.0)      { _matPtr=blosum80mtVec;    }
            else if(pcid>60.0)     { _matPtr=blosum62mt2Vec;   }
            else if(pcid>40.0)     { _matPtr=blosum45mtVec;    }
            else if(pcid>30.0)     { scaleParam.scale=0.5f; _matPtr=blosum45mtVec; }
            else if(pcid>20.0)     { scaleParam.scale=0.6f; _matPtr=blosum45mtVec; }
            else
            {
                scaleParam.scale=0.6f; _matPtr=blosum30mtVec;
            }
            _matXref=&defaultAAXref;
        }
        else if(matrixName->compare("pam")==0)
        {
            scaleParam.scale=0.75f;
            if(_negMatrix || !userParameters->getDistanceTree())
            {
                _matPtr=pam120mtVec;
            }
            else if(pcid>80.0)     { _matPtr=pam20mtVec;   }
            else if(pcid>60.0)     { _matPtr=pam60mtVec;   }
            else if(pcid>40.0)     { _matPtr=pam120mtVec;  }
            else                   { _matPtr=pam350mtVec;  }
            _matXref=&defaultAAXref;
        }
        else if(matrixName->compare("gonnet")==0)
        {
            scaleParam.scale/=2.0f;
            if(_negMatrix || !userParameters->getDistanceTree())
            {
                _matPtr=gon250mtVec;
            }
            else if(pcid>35.0)
            {
                _matPtr=gon80mtVec;
                scaleParam.scale/=2.0f;
            }
            else if(pcid>25.0)
            {
                if(minLen<100) _matPtr=gon250mtVec;
                else _matPtr=gon120mtVec;
            }
            else
            {
                if(minLen<100) _matPtr=gon350mtVec;
                else _matPtr=gon160mtVec;
            }
            _matXref=&defaultAAXref;
            scaleParam.intScale/=10;
        }
        else if(matrixName->compare("id")==0)
        {
            _matPtr=idmatVec;
            _matXref=&defaultAAXref;
        }
        else if(userSeries)
        {
            bool found=false; int j=0;
            for(int i=0;i<matSeries.nmat;i++)
            {
                if(pcid>=matSeries.mat[i].llimit && pcid<=matSeries.mat[i].ulimit)
                {
                    j=i;found=true;break;
                }
            }
            if(!found)
            {
                utilityObject->warning(
                  "Series matrix not found for pcid=%.1f.\nUsing first matrix in series.\n", pcid);
                j=0;
            }
            _matPtr=matSeries.mat[j].matptr;
            _matXref=matSeries.mat[j].AAXref;
            scaleParam.scale = (float)(0.5 +
                (pcid - matSeries.mat[j].llimit)/
                ((matSeries.mat[j].ulimit - matSeries.mat[j].llimit)*2.0));
        }
        else
        {
            // userMat
            _matPtr=&userMat;
            _matXref=&AAXref;
        }

        maxRes=getMatrix(_matPtr,_matXref,matrix,_negMatrix,(int)scaleParam.intScale);
        if(maxRes==0)
        {
            cout<<"Error: matrix "<<*matrixName<<" not found\n";
            return -1;
        }
    }

#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        ostringstream outs;
        outs<<"    Called getMatrix.\n intScale="<<scaleParam.intScale
            <<", scale="<<scaleParam.scale<<", pcid="<<pcid<<"\n\n";
        logObject->logMsg(outs.str());
    }
#endif

    matAvg=matrixAvgScore;
    return maxRes;
}

/**
 * getMatrix() – função que efetivamente carrega a submatriz
 */
int SubMatrix::getMatrix(Matrix *matptr, Xref *xref,
                         int matrix[NUMRES][NUMRES],
                         bool negFlag, int scale,
                         bool minimise)
{
    int i, j, k, ix=0;
    int ti, tj, maxRes=0;
    int av1,av2,av3;
    int minVal, maxVal;

    for(i=0;i<NUMRES;i++)
        for(j=0;j<NUMRES;j++)
            matrix[i][j]=0;

    // Preenche matrix
    ix=0;
    for(i=0; i<=userParameters->getMaxAA(); i++)
    {
        ti=(*xref)[i];
        for(j=0; j<=i; j++)
        {
            tj=(*xref)[j];
            if((ti!=-1)&&(tj!=-1))
            {
                k=(*matptr)[ix];
                if(ti==tj)
                {
                    matrix[ti][ti]=k*scale;
                    maxRes++;
                }
                else
                {
                    matrix[ti][tj]=k*scale;
                    matrix[tj][ti]=k*scale;
                }
                ix++;
            }
        }
    }

    --maxRes; // pois contava diagonal

    // av1, av2, av3
    av1=av2=av3=0;
    for(i=0;i<=userParameters->getMaxAA();i++)
    {
        for(j=0;j<=i;j++)
        {
            av1 += matrix[i][j];
            if(i==j) av2+=matrix[i][j];
            else av3+=matrix[i][j];
        }
    }

    av1 /= (maxRes*maxRes)/2; 
    av2 /= maxRes;
    av3  = (int)(av3 /(((float)(maxRes*maxRes - maxRes))/2.0));
    matrixAvgScore = -av3; 

    // encontrar minVal e maxVal
    minVal = maxVal = matrix[0][0];
    for(i=0;i<=userParameters->getMaxAA();i++)
    {
        for(j=1;j<=i;j++)
        {
            if(matrix[i][j]<minVal) minVal=matrix[i][j];
            if(matrix[i][j]>maxVal) maxVal=matrix[i][j];
        }
    }

    if(!minimise)
    {
        // se negFlag=false => tornar tudo >=0
        if(!negFlag)
        {
            if(minVal<0)
            {
                for(i=0;i<=userParameters->getMaxAA();i++)
                {
                    ti=(*xref)[i];
                    if(ti!=-1)
                    {
                        for(j=0;j<=userParameters->getMaxAA();j++)
                        {
                            tj=(*xref)[j];
                            if(tj!=-1)
                            {
                                matrix[ti][tj]-=minVal;
                            }
                        }
                    }
                }
            }
        }

        // Ajustar GAP (score=0)
        int ggScore=0, grScore=0;
        int _gapPos1=userParameters->getGapPos1();
        int _gapPos2=userParameters->getGapPos2();

        for(i=0;i<_gapPos1;i++)
        {
            matrix[i][_gapPos1]=grScore;
            matrix[_gapPos1][i]=grScore;
            matrix[i][_gapPos2]=grScore;
            matrix[_gapPos2][i]=grScore;
        }
        matrix[_gapPos1][_gapPos1]=ggScore;
        matrix[_gapPos2][_gapPos2]=ggScore;
        matrix[_gapPos2][_gapPos1]=ggScore;
        matrix[_gapPos1][_gapPos2]=ggScore;
    }
    else
    {
        // Minimizar => Saga matrix
        for(i=0;i<=userParameters->getMaxAA();i++)
        {
            for(j=0;j<=userParameters->getMaxAA();j++)
            {
                matrix[i][j] = maxVal - matrix[i][j];
            }
        }
    }

    maxRes+=2;
    return maxRes;
}

/**
 * Lê uma ou série de matrizes para user usage.
 */
bool SubMatrix::getUserMatFromFile(char *str, int alignResidueType, int alignType)
{
    checkResidueAndAlignType(alignResidueType, alignType);

    if(userParameters->getMenuFlag())
    {
        utilityObject->getStr("Enter name of the matrix file", line2);
    }
    else
    {
        line2=string(str);
    }
    if(line2.size()==0) return false;

    FILE* infile=fopen(line2.c_str(),"r");
    if(!infile)
    {
        utilityObject->error("Cannot find matrix file [%s]", line2.c_str());
        return false;
    }
    fclose(infile);
    strcpy(str, line2.c_str());

    // descobrir qual mat/xref iremos usar
    mat  = getUserMatAddress(alignResidueType, alignType);
    xref = getUserXrefAddress(alignResidueType, alignType);

    int maxRes=0;
    // se for Protein + MultipleAlign => tenta ler series:
    if((alignResidueType==Protein)&&(alignType==MultipleAlign))
    {
        maxRes=readMatrixSeries(str, userMat, AAXref);
    }
    else
    {
        maxRes=readUserMatrix(str, *mat, *xref);
    }
    return (maxRes>0);
}

bool SubMatrix::getUserMatSeriesFromFile(char *str)
{
    if(userParameters->getMenuFlag())
    {
        utilityObject->getStr("Enter name of the matrix file", line2);
    }
    else
    {
        line2=string(str);
    }
    if(line2.size()==0) return false;

    FILE*infile=fopen(line2.c_str(),"r");
    if(!infile)
    {
        utilityObject->error("Cannot find matrix file [%s]", line2.c_str());
        return false;
    }
    fclose(infile);
    strcpy(str,line2.c_str());

    int maxRes=readMatrixSeries(str, userMat, AAXref);
    return (maxRes>0);
}

/**
 * Lê várias matrizes (CLUSTAL_SERIES).
 */
int SubMatrix::readMatrixSeries(const char *fileName, Matrix &userMat, Xref &xref)
{
    FILE* fd=fopen(fileName,"r");
    if(!fd)
    {
        utilityObject->error("cannot open %s",fileName);
        return 0;
    }
    char inline1[1024];
    userSeries=false;
    while(fgets(inline1,1024,fd))
    {
        if(commentline(inline1)) continue;
        if(utilityObject->lineType(inline1,"CLUSTAL_SERIES"))
            userSeries=true;
        else
            userSeries=false;
        break;
    }

    if(!userSeries)
    {
        fclose(fd);
        int r=readUserMatrix(fileName, userMat, xref);
        return r;
    }

    // se é series
    int nmat=0; matSeries.nmat=0;
    rewind(fd);
    while(fgets(inline1,1024,fd))
    {
        if(commentline(inline1)) continue;
        if(utilityObject->lineType(inline1,"MATRIX"))
        {
            int llimit, ulimit;
            char mat_fileName[FILENAMELEN];
            if(sscanf(inline1+6,"%d %d %s",&llimit,&ulimit,mat_fileName)!=3)
            {
                utilityObject->error("Bad format in file %s\n",fileName);
                fclose(fd);return 0;
            }
            if(llimit<0||llimit>100||ulimit<0||ulimit>100)
            {
                utilityObject->error("Bad format in file %s\n",fileName);
                fclose(fd);return 0;
            }
            if(ulimit<=llimit)
            {
                utilityObject->error("in file %s: lower>upper (%d-%d)\n",fileName,llimit,ulimit);
                fclose(fd);return 0;
            }
            int rr=readUserMatrix(mat_fileName, userMatSeries[nmat], AAXrefseries[nmat]);
            if(rr<=0)
            {
                utilityObject->error("Bad format in matrix file %s\n",mat_fileName);
                fclose(fd);return 0;
            }
            matSeries.mat[nmat].llimit=llimit;
            matSeries.mat[nmat].ulimit=ulimit;
            matSeries.mat[nmat].matptr=&userMatSeries[nmat];
            matSeries.mat[nmat].AAXref=&AAXrefseries[nmat];
            nmat++;
            if(nmat>=MAXMAT)
            {
                cout<<"Matrix series file has more entries than allowed (MAXMAT="<<MAXMAT
                    <<"). Using first "<<MAXMAT<<" only.\n";
                break;
            }
        }
    }
    fclose(fd);
    matSeries.nmat=nmat;
    return (nmat>0 ? 20 : 0); // se nmat>0 devolve algo>0
}

/**
 * Lê uma única matriz do arquivo.
 */
int SubMatrix::readUserMatrix(const char *fileName, Matrix &userMat, Xref &xref)
{
    if(fileName[0]=='\0')
    {
        utilityObject->error("comparison matrix not specified");
        return 0;
    }
    FILE* fd=fopen(fileName,"r");
    if(!fd)
    {
        utilityObject->error("cannot open %s",fileName);
        return 0;
    }

    char inline1[1024];
    char codes[NUMRES];
    int k=0;

    while(fgets(inline1,1024,fd))
    {
        if(commentline(inline1)) continue;
        if(utilityObject->lineType(inline1,"CLUSTAL_SERIES"))
        {
            utilityObject->error("in %s - single matrix expected.",fileName);
            fclose(fd);
            return 0;
        }
        // ler chars
        k=0;
        for(int idx=0; idx<(int)strlen(inline1); idx++)
        {
            if(isalpha((unsigned char)inline1[idx]))
            {
                codes[k++]=inline1[idx];
                if(k>NUMRES)
                {
                    utilityObject->error("too many entries in matrix %s",fileName);
                    fclose(fd);
                    return 0;
                }
            }
        }
        codes[k]='\0';
        break;
    }
    if(k==0)
    {
        utilityObject->error("wrong format in matrix %s",fileName);
        fclose(fd);
        return 0;
    }

    // crossref
    for(int i=0;i<NUMRES;i++)
        xref[i]=-1;

    int maxRes=0;
    for(int i=0;i<k;i++)
    {
        char c1=codes[i];
        bool found=false;
        for(int j=0; ; j++)
        {
            char c2=userParameters->getAminoAcidCode(j);
            if(!c2) break; // fim
            if(c1==c2)
            {
                xref[i]=j;
                maxRes++;
                found=true;
                break;
            }
        }
        if(!found && c1!='*')
        {
            utilityObject->warning("residue %c in matrix %s not recognised", c1, fileName);
        }
    }

    int ix=0, ix1=0;
    while(fgets(inline1,1024,fd))
    {
        if(inline1[0]=='\n') continue;
        if(inline1[0]=='#' || inline1[0]=='!')
            break;
        char* args[NUMRES+4];
        int numargs = getArgs(inline1,args,k+1);
        if(numargs<maxRes)
        {
            utilityObject->error("wrong format in matrix %s",fileName);
            fclose(fd);
            return 0;
        }
        int farg=0;
        if(isalpha((unsigned char)args[0][0])) farg=1;

        float scale=1.0f;
        for(unsigned t=0;t<strlen(args[farg]);t++)
        {
            if(args[farg][t]=='.') { scale=10.0f; break; }
        }
        for(int i=0;i<=ix;i++)
        {
            if(xref[i]!=-1)
            {
                float ff=atof(args[i+farg]);
                userMat[ix1++]=(short)(ff*scale);
            }
        }
        ix++;
    }

    if(ix!=(k+1))
    {
        utilityObject->error("wrong format in matrix %s",fileName);
        fclose(fd);
        return 0;
    }

    userMat.resize(ix1+1);
    fclose(fd);

    maxRes+=2;
    return maxRes;
}

/**
 * compareMatrices() - debug
 */
void SubMatrix::compareMatrices(int mat1[NUMRES][NUMRES],
                                int mat2[NUMRES][NUMRES])
{
    int same=1;
    for(int r=0;r<NUMRES;r++)
    {
        for(int c=0;c<NUMRES;c++)
        {
            if(mat1[r][c]!=mat2[r][c])
            {
                same=0;
                cout<<"Row="<<r<<",Col="<<c<<" difere.\n";
                break;
            }
        }
        if(!same) break;
    }
    if(!same) cout<<"As matrizes são diferentes!\n";
    else      cout<<"As matrizes são iguais!\n";
}

/**
 * printGetMatrixResults() - debug
 */
void SubMatrix::printGetMatrixResults(int mat[NUMRES][NUMRES])
{
    ofstream outfile("getmatrix.out");
    if(!outfile)
    {
        cerr<<"falha ao abrir getmatrix.out!\n";
        return;
    }
    for(int row=0; row<NUMRES; row++)
    {
        for(int col=0; col<NUMRES; col++)
        {
            if(mat[row][col]>9 || mat[row][col]<0)
                outfile<<" "<<mat[row][col]<<",";
            else
                outfile<<"  "<<mat[row][col]<<",";
        }
        outfile<<"\n";
    }
}

/**
 * getAlnScoreMatrix(): Se quiser BLOSUM45 p/ scoring final
 */
int SubMatrix::getAlnScoreMatrix(int matrix[NUMRES][NUMRES])
{
    int _maxNumRes = getMatrix(blosum45mtVec, &defaultAAXref,
                               matrix, true, 100);
    return _maxNumRes;
}

/**
 * getQTMatrixForHistogram() – ClustalQt
 */
void SubMatrix::getQTMatrixForHistogram(int matrix[NUMRES][NUMRES])
{
    Matrix* matPtrLocal;
    Xref*   xrefLocal;
    int maxRes=0;

    if(userParameters->getDNAFlag())
    {
        if(QTDNAHistMatNum==DNAUSERDEFINED)
        {
            matPtrLocal=&QTscoreUserDNAMatrix;
            xrefLocal=&QTscoreDNAXref;
        }
        else if(QTDNAHistMatNum==DNACLUSTALW)
        {
            matPtrLocal=clustalvdnamtVec;
            xrefLocal=&defaultDNAXref;
        }
        else
        {
            matPtrLocal=swgapdnamtVec;
            xrefLocal=&defaultDNAXref;
        }
    }
    else
    {
        if(QTAAHistMatNum==AAHISTIDENTITY)
        {
            matPtrLocal=idmatVec;
            xrefLocal=&defaultAAXref;
        }
        else if(QTAAHistMatNum==AAHISTGONNETPAM80)
        {
            matPtrLocal=gon80mtVec;
            xrefLocal=&defaultAAXref;
        }
        else if(QTAAHistMatNum==AAHISTGONNETPAM120)
        {
            matPtrLocal=gon120mtVec;
            xrefLocal=&defaultAAXref;
        }
        else if(QTAAHistMatNum==AAHISTUSER)
        {
            matPtrLocal=&QTscoreUserMatrix;
            xrefLocal=&QTscoreXref;
        }
        else if(QTAAHistMatNum==AAHISTGONNETPAM350)
        {
            matPtrLocal=gon350mtVec;
            xrefLocal=&defaultAAXref;
        }
        else
        {
            // default gon250
            matPtrLocal=gon250mtVec;
            xrefLocal=&defaultAAXref;
        }
    }
    maxRes=getMatrix(matPtrLocal,xrefLocal,matrix,false,100);
}

/**
 * getQTMatrixForLowScoreSeg() – ClustalQt
 */
void SubMatrix::getQTMatrixForLowScoreSeg(int matrix[NUMRES][NUMRES])
{
    Matrix* matPtrLocal;
    Xref*   xrefLocal;
    int maxRes=0;
    int _maxAA=userParameters->getMaxAA();
    int maxv=0;

    if(userParameters->getDNAFlag())
    {
        if(QTsegmentDNAMatNum==DNAUSERDEFINED)
        {
            matPtrLocal=&QTsegmentDNAMatrix;
            xrefLocal=&QTsegmentDNAXref;
        }
        else if(QTsegmentDNAMatNum==DNACLUSTALW)
        {
            matPtrLocal=clustalvdnamtVec;
            xrefLocal=&defaultDNAXref;
        }
        else
        {
            matPtrLocal=swgapdnamtVec;
            xrefLocal=&defaultDNAXref;
        }
        maxRes=getMatrix(matPtrLocal,xrefLocal,matrix,false,100);

        for(int i=0;i<=_maxAA;i++)
            for(int j=0;j<=_maxAA;j++)
                if(matrix[i][j]>maxv) maxv=matrix[i][j];

        int offset=(int)( (float)maxv * userParameters->getQTlowScoreDNAMarkingScale()/20.0f );
        for(int i=0;i<=_maxAA;i++)
            for(int j=0;j<=_maxAA;j++)
                matrix[i][j]-=offset;
    }
    else
    {
        if(QTsegmentAAMatNum==QTAASEGGONNETPAM80)
        {
            matPtrLocal=gon80mtVec;
            xrefLocal=&defaultAAXref;
        }
        else if(QTsegmentAAMatNum==QTAASEGGONNETPAM120)
        {
            matPtrLocal=gon120mtVec;
            xrefLocal=&defaultAAXref;
        }
        else if(QTsegmentAAMatNum==QTAASEGUSER)
        {
            matPtrLocal=&QTsegmentAAMatrix;
            xrefLocal=&QTsegmentAAXref;
        }
        else if(QTsegmentAAMatNum==QTAASEGGONNETPAM350)
        {
            matPtrLocal=gon350mtVec;
            xrefLocal=&defaultAAXref;
        }
        else
        {
            matPtrLocal=gon250mtVec;
            xrefLocal=&defaultAAXref;
        }
        // get negative matrix
        maxRes=getMatrix(matPtrLocal,xrefLocal,matrix,true,100);
    }
}

bool SubMatrix::getQTLowScoreMatFromFile(char* fileName, bool dna)
{
    FILE*infile=fopen(fileName,"r");
    if(!infile)
    {
        utilityObject->error("Cannot find matrix file [%s]",fileName);
        return false;
    }
    fclose(infile);

    if(dna)
    {
        int r=readUserMatrix(fileName,QTsegmentDNAMatrix,QTsegmentDNAXref);
        return (r>0);
    }
    else
    {
        int r=readUserMatrix(fileName,QTsegmentAAMatrix,QTsegmentAAXref);
        return (r>0);
    }
}

bool SubMatrix::getAAScoreMatFromFile(char* str)
{
    FILE*infile=fopen(str,"r");
    if(!infile)
    {
        utilityObject->error("Cannot find matrix file [%s]",str);
        return false;
    }
    fclose(infile);
    int r=readUserMatrix(str,QTscoreUserMatrix,QTscoreXref);
    return (r>0);
}

bool SubMatrix::getDNAScoreMatFromFile(char* str)
{
    FILE*infile=fopen(str,"r");
    if(!infile)
    {
        utilityObject->error("Cannot find matrix file [%s]",str);
        return false;
    }
    fclose(infile);
    int r=readUserMatrix(str,QTscoreUserDNAMatrix,QTscoreDNAXref);
    return (r>0);
}

int SubMatrix::getMatrixNum()      { return matrixNum; }
int SubMatrix::getDNAMatrixNum()   { return DNAMatrixNum; }
int SubMatrix::getPWMatrixNum()    { return pwMatrixNum; }
int SubMatrix::getPWDNAMatrixNum() { return pwDNAMatrixNum; }

/**
 * setCurrentNameAndNum() - define qual matrix ou series usar
 */
void SubMatrix::setCurrentNameAndNum(string _matrixName, int _matrixNum,
                                     int alignResidueType,int alignType)
{
    checkResidueAndAlignType(alignResidueType, alignType);

    if((alignResidueType==Protein)&&(alignType==Pairwise))
    {
        pwMatrixNum=_matrixNum;
        delete pwMatrixName; // se já existia
        pwMatrixName=new string(_matrixName);
    }
    else if((alignResidueType==Protein)&&(alignType==MultipleAlign))
    {
        matrixNum=_matrixNum;
        delete matrixName;
        matrixName=new string(_matrixName);
    }
    else if((alignResidueType==DNA)&&(alignType==Pairwise))
    {
        pwDNAMatrixNum=_matrixNum;
        delete pwDNAMatrixName;
        pwDNAMatrixName=new string(_matrixName);
    }
    else if((alignResidueType==DNA)&&(alignType==MultipleAlign))
    {
        DNAMatrixNum=_matrixNum;
        delete DNAMatrixName;
        DNAMatrixName=new string(_matrixName);
    }
#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        ostringstream outs;
        outs<<"setCurrentNameAndNum: (residueType="<<alignResidueType
            <<",alignType="<<alignType<<") => "<<_matrixName<<"\n";
        logObject->logMsg(outs.str());
    }
#endif
}

/**
 * getMatrixNumForMenu()
 */
int SubMatrix::getMatrixNumForMenu(int alignResidueType,int alignType)
{
    checkResidueAndAlignType(alignResidueType,alignType);
    if((alignResidueType==Protein)&&(alignType==Pairwise))
        return pwMatrixNum;
    else if((alignResidueType==Protein)&&(alignType==MultipleAlign))
        return matrixNum;
    else if((alignResidueType==DNA)&&(alignType==Pairwise))
        return pwDNAMatrixNum;
    else if((alignResidueType==DNA)&&(alignType==MultipleAlign))
        return DNAMatrixNum;
    return -100; // fallback
}

/**
 * commentline() - auxiliar
 */
bool SubMatrix::commentline(char* line)
{
    if(line[0]=='#') return true;
    for(int i=0; line[i] && line[i]!= '\n'; i++)
    {
        if(!isspace((unsigned char)line[i])) return false;
    }
    return true;
}

/**
 * printInFormat() — recebe “const char* name”
 */
void SubMatrix::printInFormat(vector<short>& temp, const char* name)
{
    char nameFile[60];
    sprintf(nameFile, "%s.out", name);

    ofstream outfile(nameFile);
    if(!outfile)
    {
        cerr<<"falha ao abrir "<<nameFile<<"\n";
        return;
    }
    outfile<<"short "<<name<<"[]{\n";
    int numOnLine=1;
    int soFar=0;
    for(int i=0; i<(int)temp.size(); i++)
    {
        if(soFar==numOnLine)
        {
            outfile<<"\n";
            soFar=0;
            numOnLine++;
        }
        if(temp[i]>9 || temp[i]<0) outfile<<" "<<temp[i]<<",";
        else                      outfile<<"  "<<temp[i]<<",";
        soFar++;
        if(i+1==(int)temp.size()-1)
        {
            if(temp[i+1]>9||temp[i+1]<0)
                outfile<<" "<<temp[i+1]<<"};\n";
            else
                outfile<<"  "<<temp[i+1]<<"};\n";
            break;
        }
    }
    // Em paralelo, “temp.out”
    ofstream outfile2("temp.out");
    for(size_t i=0; i<temp.size(); i++)
    {
        outfile2<<temp[i]<<" ";
    }
}

/**
 * printVectorToFile() — também usa “const char* name”
 */
void SubMatrix::printVectorToFile(vector<short>& temp, const char* name)
{
    char nameFile[60];
    sprintf(nameFile, "%s.out", name);

    ofstream outfile(nameFile);
    if(!outfile)
    {
        cerr<<"falha ao abrir "<<nameFile<<"\n";
        return;
    }
    for(int i=0; i<(int)temp.size(); i++)
    {
        if(temp[i]>9||temp[i]<0) outfile<<" "<<temp[i]<<",";
        else                    outfile<<"  "<<temp[i]<<",";
    }
    outfile.close();
}

/**
 * getUserMatAddress() — retorna ref p/ userMat/pwUserMat etc
 */
Matrix* SubMatrix::getUserMatAddress(int alignResidueType, int alignType)
{
    if((alignResidueType==Protein)&&(alignType==Pairwise))
        return &pwUserMat;
    else if((alignResidueType==Protein)&&(alignType==MultipleAlign))
        return &userMat;
    else if((alignResidueType==DNA)&&(alignType==Pairwise))
        return &pwUserDNAMat;
    else
        return &userDNAMat; // DNA+Multiple
}

/**
 * getUserXrefAddress()
 */
Xref* SubMatrix::getUserXrefAddress(int alignResidueType, int alignType)
{
    if((alignResidueType==Protein)&&(alignType==Pairwise))
        return &pwAAXref;
    else if((alignResidueType==Protein)&&(alignType==MultipleAlign))
        return &AAXref;
    else if((alignResidueType==DNA)&&(alignType==Pairwise))
        return &pwDNAXref;
    else
        return &DNAXref;
}

/**
 * checkResidueAndAlignType() - se valores incorretos, erro e sai
 */
void SubMatrix::checkResidueAndAlignType(int alignResidueType, int alignType)
{
    // 0=Protein,1=DNA; 0=Pairwise,1=MultipleAlign
    if(((alignResidueType!=0)&&(alignResidueType!=1))
       || ((alignType!=0)&&(alignType!=1)))
    {
        InvalidCombination ex(alignResidueType,alignType);
        ex.whatHappened();
        exit(1);
    }
}

/**
 * Rotina de teste
 */
void SubMatrix::tempInterface(int alignResidueType,int alignType)
{
    int matrix[NUMRES][NUMRES];
    PairScaleValues temp;
    userParameters->setDNAFlag(true);
    PrfScaleValues scaleParam;
    char userFile[FILENAMELEN+1];
    strcpy(userFile,"mat1");
    userParameters->setMenuFlag(false);

    getUserMatFromFile(userFile,DNA,Pairwise);
    setCurrentNameAndNum(userFile,4,DNA,Pairwise);
    setCurrentNameAndNum("gonnet",4,Protein,Pairwise);
}

} // end namespace clustalw
