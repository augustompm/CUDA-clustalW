/**
 * Author: Mark Larkin
 *
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.
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
#include "SubMatrix.h"
#include "matrices.h"
#include "../general/InvalidCombination.cpp"
#include "../general/debuglogObject.h"

namespace clustalw
{
using namespace std;

/**
 * Construtor SubMatrix com checagens mínimas para variáveis.
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
        blosum30mtVec    = new Matrix(blosum30mt,    blosum30mt    + sizenAAMatrix);
        blosum40mtVec    = new Matrix(blosum40mt,    blosum40mt    + sizenAAMatrix);
        blosum45mtVec    = new Matrix(blosum45mt,    blosum45mt    + sizenAAMatrix);
        blosum62mt2Vec   = new Matrix(blosum62mt2,   blosum62mt2   + sizenAAMatrix);
        blosum80mtVec    = new Matrix(blosum80mt,    blosum80mt    + sizenAAMatrix);
        pam20mtVec       = new Matrix(pam20mt,       pam20mt       + sizenAAMatrix);
        pam60mtVec       = new Matrix(pam60mt,       pam60mt       + sizenAAMatrix);
        pam120mtVec      = new Matrix(pam120mt,      pam120mt      + sizenAAMatrix);
        pam350mtVec      = new Matrix(pam350mt,      pam350mt      + sizenAAMatrix);
        idmatVec         = new Matrix(idmat,         idmat         + sizenAAMatrix);
        gon40mtVec       = new Matrix(gon40mt,       gon40mt       + sizenAAMatrix);
        gon80mtVec       = new Matrix(gon80mt,       gon80mt       + sizenAAMatrix);
        gon120mtVec      = new Matrix(gon120mt,      gon120mt      + sizenAAMatrix);
        gon160mtVec      = new Matrix(gon160mt,      gon160mt      + sizenAAMatrix);
        gon250mtVec      = new Matrix(gon250mt,      gon250mt      + sizenAAMatrix);
        gon350mtVec      = new Matrix(gon350mt,      gon350mt      + sizenAAMatrix);

        clustalvdnamtVec = new Matrix(clustalvdnamt, clustalvdnamt + sizeDNAMatrix);
        swgapdnamtVec    = new Matrix(swgapdnamt,    swgapdnamt    + sizeDNAMatrix);

        userMat.resize(NUMRES * NUMRES);
        pwUserMat.resize(NUMRES * NUMRES);
        userDNAMat.resize(NUMRES * NUMRES);
        pwUserDNAMat.resize(NUMRES * NUMRES);
        QTscoreUserMatrix.resize(NUMRES * NUMRES);
        QTscoreUserDNAMatrix.resize(NUMRES * NUMRES);
        QTsegmentDNAMatrix.resize(NUMRES * NUMRES);
        QTsegmentAAMatrix.resize(NUMRES * NUMRES);

        userMatSeries.resize(MAXMAT);
        {
            auto firstM = userMatSeries.begin();
            auto lastM  = userMatSeries.end();
            while (firstM != lastM)
            {
                firstM->resize(NUMRES * NUMRES);
                ++firstM;
            }
        }

        AAXrefseries.resize(MAXMAT);
        {
            auto firstX = AAXrefseries.begin();
            auto lastX  = AAXrefseries.end();
            while (firstX != lastX)
            {
                firstX->resize(NUMRES + 1);
                ++firstX;
            }
        }

        matrixNum       = 3;
        matrixName      = new string("gonnet");
        DNAMatrixNum    = 1;
        DNAMatrixName   = new string("iub");
        pwMatrixNum     = 3;
        pwMatrixName    = new string("gonnet");
        pwDNAMatrixNum  = 1;
        pwDNAMatrixName = new string("iub");
    }
    catch(const exception &ex)
    {
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        exit(1);
    }
}

void SubMatrix::setValuesToDefault()
{
    matrixAvgScore      = 0;
    QTDNAHistMatNum     = DNAIUB;
    QTAAHistMatNum      = AAHISTGONNETPAM250;
    QTsegmentDNAMatNum  = DNAIUB;
    QTsegmentAAMatNum   = QTAASEGGONNETPAM250;
    userSeries          = false;
    matrixNum           = 3;
    DNAMatrixNum        = 1;
    pwMatrixNum         = 3;
    pwDNAMatrixNum      = 1;
}

SubMatrix::~SubMatrix()
{
    // Liberação da memória alocada dinamicamente
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

    delete matrixName;
    delete DNAMatrixName;
    delete pwMatrixName;
    delete pwDNAMatrixName;
}

/**
 * Configura cross references para aminoácidos e DNA.
 */
void SubMatrix::setUpCrossReferences()
{
    char c1, c2;
    short i, j, maxRes;
    defaultAAXref.resize(NUMRES + 1);
    defaultDNAXref.resize(NUMRES + 1);

    string aminoAcidOrder   = "ABCDEFGHIKLMNPQRSTVWXYZ";
    string nucleicAcidOrder = "ABCDGHKMNRSTUVWXY";

    DNAXref.resize(NUMRES + 1);
    AAXref.resize(NUMRES + 1);
    pwAAXref.resize(NUMRES + 1);
    pwDNAXref.resize(NUMRES + 1);
    QTscoreXref.resize(NUMRES + 1);
    QTscoreDNAXref.resize(NUMRES + 1);
    QTsegmentDNAXref.resize(NUMRES + 1);
    QTsegmentAAXref.resize(NUMRES + 1);

    for (i = 0; i < NUMRES; i++)
    {
        defaultAAXref[i]  = -1;
        defaultDNAXref[i] = -1;
    }

    maxRes = 0;
    for (i = 0; (c1 = aminoAcidOrder[i]); i++)
    {
        for (j = 0; (c2 = userParameters->getAminoAcidCode(j)); j++)
        {
            if (c1 == c2)
            {
                defaultAAXref[i] = j;
                maxRes++;
                break;
            }
        }
        if ((defaultAAXref[i] == -1) && (aminoAcidOrder[i] != '*'))
        {
            utilityObject->error("residue %c in matrices.h is not recognised",
                                 aminoAcidOrder[i]);
        }
    }

    maxRes = 0;
    for (i = 0; (c1 = nucleicAcidOrder[i]); i++)
    {
        for (j = 0; (c2 = userParameters->getAminoAcidCode(j)); j++)
        {
            if (c1 == c2)
            {
                defaultDNAXref[i] = j;
                maxRes++;
                break;
            }
        }
        if ((defaultDNAXref[i] == -1) && (nucleicAcidOrder[i] != '*'))
        {
            utilityObject->error("nucleic acid %c in matrices.h is not recognised",
                                 nucleicAcidOrder[i]);
        }
    }
}

int SubMatrix::getPairwiseMatrix(int matrix[NUMRES][NUMRES],
                                 PairScaleValues& scale,
                                 int& matAvg)
{
    if (!this)
    {
        throw std::runtime_error("SubMatrix::getPairwiseMatrix() - 'this' é inválido (nulo).");
    }

    int _maxRes;
    Matrix* _matPtr;
    Xref*   _matXref;

#ifdef OS_MAC
    scale.intScale = 10;
#else
    scale.intScale = 100;
#endif

    scale.gapOpenScale   = 1.0;
    scale.gapExtendScale = 1.0;

#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        logObject->logMsg("In the function getPairwiseMatrix:\n");
    }
#endif

    if (userParameters->getDNAFlag())
    {
#if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            string msg = "    (DNA AND Pairwise) " + *pwDNAMatrixName + "\n";
            logObject->logMsg(msg);
        }
#endif
        if (pwDNAMatrixName->compare("iub") == 0)
        {
            _matPtr  = swgapdnamtVec;
            _matXref = &defaultDNAXref;
        }
        else if (pwDNAMatrixName->compare("clustalw") == 0)
        {
            _matPtr  = clustalvdnamtVec;
            _matXref = &defaultDNAXref;
            scale.gapOpenScale   = 0.6667f;
            scale.gapExtendScale = 0.751f;
        }
        else
        {
            _matPtr  = &pwUserDNAMat;
            _matXref = &pwDNAXref;
        }
        _maxRes = getMatrix(_matPtr, _matXref, matrix, true, scale.intScale);

        if (_maxRes == 0)
        {
            return -1;
        }
        float _transitionWeight = userParameters->getTransitionWeight();
        // Ajustes no matrix de transição
        matrix[0][4]    = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[4][0]    = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[2][11]   = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[11][2]   = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[2][12]   = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[12][2]   = static_cast<int>(_transitionWeight * matrix[0][0]);
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
        if (pwMatrixName->compare("blosum") == 0)
        {
            _matPtr  = blosum30mtVec;
            _matXref = &defaultAAXref;
        }
        else if (pwMatrixName->compare("pam") == 0)
        {
            _matPtr  = pam350mtVec;
            _matXref = &defaultAAXref;
        }
        else if (pwMatrixName->compare("gonnet") == 0)
        {
            _matPtr  = gon250mtVec;
            _matXref = &defaultAAXref;
            scale.intScale /= 10;
        }
        else if (pwMatrixName->compare("id") == 0)
        {
            _matPtr  = idmatVec;
            _matXref = &defaultAAXref;
        }
        else
        {
            _matPtr  = &pwUserMat;
            _matXref = &pwAAXref;
        }

        _maxRes = getMatrix(_matPtr, _matXref, matrix, true, scale.intScale);
        if (_maxRes == 0)
        {
            return -1;
        }
    }

#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        ostringstream outs;
        outs << "    Called getMatrix with ?  intScale=" << scale.intScale
             << ", gapOpenScale=" << scale.gapOpenScale
             << ", gapExtendScale=" << scale.gapExtendScale
             << "\n\n";
        logObject->logMsg(outs.str());
    }
#endif

    matAvg = matrixAvgScore;
    return _maxRes;
}

int SubMatrix::getProfileAlignMatrix(int matrix[NUMRES][NUMRES],
                                     double pcid, int minLen,
                                     PrfScaleValues& scaleParam,
                                     int& matAvg)
{
    if (!this)
    {
        throw std::runtime_error("SubMatrix::getProfileAlignMatrix() - 'this' é inválido (nulo).");
    }

    bool  _negMatrix = userParameters->getUseNegMatrix();
    int   _maxRes    = 0;
    scaleParam.intScale = 100;
    scaleParam.scale    = 0.75f;
    Matrix* _matPtr;
    Xref*   _matXref;
    
#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        logObject->logMsg("In the function getProfileAlignMatrix:\n");
    }
#endif

    if (userParameters->getDNAFlag())
    {
#if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            ostringstream outs;
            outs << "    (DNA AND Multiple align) " << DNAMatrixName->c_str() << "\n";
            logObject->logMsg(outs.str());
        }
#endif
        scaleParam.scale = 1.0f;
        if (DNAMatrixName->compare("iub") == 0)
        {
            _matPtr  = swgapdnamtVec;
            _matXref = &defaultDNAXref;
        }
        else if (DNAMatrixName->compare("clustalw") == 0)
        {
            _matPtr  = clustalvdnamtVec;
            _matXref = &defaultDNAXref;
            scaleParam.scale = 0.66f;
        }
        else
        {
            _matPtr  = &userDNAMat;
            _matXref = &DNAXref;
        }

        _maxRes = getMatrix(_matPtr, _matXref, matrix, _negMatrix,
                            scaleParam.intScale);
        if (_maxRes == 0)
        {
            return -1;
        }
        float _transitionWeight = userParameters->getTransitionWeight();
        matrix[(*_matXref)[0]][(*_matXref)[4]]  = static_cast<int>(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[4]][(*_matXref)[0]]  = static_cast<int>(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[2]][(*_matXref)[11]] = static_cast<int>(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[11]][(*_matXref)[2]] = static_cast<int>(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[2]][(*_matXref)[12]] = static_cast<int>(_transitionWeight*matrix[0][0]);
        matrix[(*_matXref)[12]][(*_matXref)[2]] = static_cast<int>(_transitionWeight*matrix[0][0]);
    }
    else
    {
#if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            ostringstream outs;
            outs << "    (Protein AND Multiple align) " << matrixName->c_str() << "\n";
            logObject->logMsg(outs.str());
        }
#endif
        if (matrixName->compare("blosum") == 0)
        {
            scaleParam.scale = 0.75f;
            if (_negMatrix || userParameters->getDistanceTree() == false)
            {
                _matPtr = blosum40mtVec;
            }
            else if (pcid > 80.0)
            {
                _matPtr = blosum80mtVec;
            }
            else if (pcid > 60.0)
            {
                _matPtr = blosum62mt2Vec;
            }
            else if (pcid > 40.0)
            {
                _matPtr = blosum45mtVec;
            }
            else if (pcid > 30.0)
            {
                scaleParam.scale = 0.5f;
                _matPtr = blosum45mtVec;
            }
            else if (pcid > 20.0)
            {
                scaleParam.scale = 0.6f;
                _matPtr = blosum45mtVec;
            }
            else
            {
                scaleParam.scale = 0.6f;
                _matPtr = blosum30mtVec;
            }
            _matXref = &defaultAAXref;
        }
        else if (matrixName->compare("pam") == 0)
        {
            scaleParam.scale = 0.75f;
            if (_negMatrix || userParameters->getDistanceTree() == false)
            {
                _matPtr = pam120mtVec;
            }
            else if (pcid > 80.0)
            {
                _matPtr = pam20mtVec;
            }
            else if (pcid > 60.0)
            {
                _matPtr = pam60mtVec;
            }
            else if (pcid > 40.0)
            {
                _matPtr = pam120mtVec;
            }
            else
            {
                _matPtr = pam350mtVec;
            }
            _matXref = &defaultAAXref;
        }
        else if (matrixName->compare("gonnet") == 0)
        {
            scaleParam.scale /= 2.0f;
            if (_negMatrix || userParameters->getDistanceTree() == false)
            {
                _matPtr = gon250mtVec;
            }
            else if (pcid > 35.0)
            {
                _matPtr = gon80mtVec;
                scaleParam.scale /= 2.0f;
            }
            else if (pcid > 25.0)
            {
                if (minLen < 100)
                {
                    _matPtr = gon250mtVec;
                }
                else
                {
                    _matPtr = gon120mtVec;
                }
            }
            else
            {
                if (minLen < 100)
                {
                    _matPtr = gon350mtVec;
                }
                else
                {
                    _matPtr = gon160mtVec;
                }
            }
            _matXref = &defaultAAXref;
            scaleParam.intScale /= 10;
        }
        else if (matrixName->compare("id") == 0)
        {
            _matPtr  = idmatVec;
            _matXref = &defaultAAXref;
        }
        else if (userSeries)
        {
            _matPtr = nullptr;
            bool found = false;
            int j = 0;
            for (int i = 0; i < matSeries.nmat; i++)
            {
                if (pcid >= matSeries.mat[i].llimit && pcid <= matSeries.mat[i].ulimit)
                {
                    j = i;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                utilityObject->warning(
                    "\nSeries matrix not found for sequence percent identity = %d.\n"
                    "(Using first matrix in series as a default.)\n"
                    "This alignment may not be optimal!\n"
                    "SUGGESTION: Check your matrix series input file and try again.",
                    (int)pcid);
                j = 0;
            }
            _matPtr  = matSeries.mat[j].matptr;
            _matXref = matSeries.mat[j].AAXref;
            scaleParam.scale = 0.5f + (pcid - matSeries.mat[j].llimit) /
                               ((matSeries.mat[j].ulimit - matSeries.mat[j].llimit) * 2.0f);
        }
        else
        {
            _matPtr  = &userMat;
            _matXref = &AAXref;
        }

        _maxRes = getMatrix(_matPtr, _matXref, matrix, _negMatrix,
                            static_cast<int>(scaleParam.intScale));
        if (_maxRes == 0)
        {
            cout << "Error: matrix " << matrixName << " not found\n";
            return -1;
        }
    }

#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        ostringstream outs;
        outs << "    Called getMatrix with ? (ProfileAlign)  intScale="
             << scaleParam.intScale << ", scale=" << scaleParam.scale
             << ", pcid=" << pcid << "\n\n";
        logObject->logMsg(outs.str());
    }
#endif

    matAvg = matrixAvgScore;
    return _maxRes;
}

/**
 * Novo getMatrix com validação extra e sem quebrar funcionalidade existente.
 */
int SubMatrix::getMatrix(Matrix* matptr,
                         Xref* xref,
                         int matrix[NUMRES][NUMRES],
                         bool negFlag,
                         int scale,
                         bool minimise)
{
    // Checagens mínimas:
    if (!this)
    {
        throw std::runtime_error("SubMatrix::getMatrix() - 'this' é inválido (nulo).");
    }
    if (!matptr)
    {
        throw std::runtime_error("SubMatrix::getMatrix() - matptr é nulo.");
    }
    if (!xref)
    {
        throw std::runtime_error("SubMatrix::getMatrix() - xref é nulo.");
    }

    int ggScore = 0;
    int grScore = 0;
    int i, j, k, ix = 0;
    int ti, tj;
    int maxRes;
    int av1, av2, av3, minVal, maxVal;

    // Inicializa a matriz
    for (i = 0; i < NUMRES; i++)
    {
        for (j = 0; j < NUMRES; j++)
        {
            matrix[i][j] = 0;
        }
    }

    ix = 0;
    maxRes = 0;

    // Preenche
    for (i = 0; i <= userParameters->getMaxAA(); i++)
    {
        ti = (*xref)[i];
        for (j = 0; j <= i; j++)
        {
            tj = (*xref)[j];
            if (ti != -1 && tj != -1)
            {
                if (ix >= static_cast<int>(matptr->size()))
                {
                    // Evita acesso fora do vetor:
                    break;
                }
                k = (*matptr)[ix];
                if (ti == tj)
                {
                    matrix[ti][ti] = k * scale;
                    maxRes++;
                }
                else
                {
                    matrix[ti][tj] = k * scale;
                    matrix[tj][ti] = k * scale;
                }
                ix++;
            }
        }
    }
    --maxRes;

    av1 = av2 = av3 = 0;
    for (i = 0; i <= userParameters->getMaxAA(); i++)
    {
        for (j = 0; j <= i; j++)
        {
            av1 += matrix[i][j];
            if (i == j)
            {
                av2 += matrix[i][j];
            }
            else
            {
                av3 += matrix[i][j];
            }
        }
    }

    if (maxRes < 2)
    {
        maxRes = 2;
    }

    av1 /= (maxRes * maxRes) / 2;
    av2 /= maxRes;
    av3 = static_cast<int>(
        av3 / ( ( static_cast<float>( maxRes * maxRes - maxRes ) ) / 2.0f )
    );
    matrixAvgScore = -av3;

    // Achar min e max
    minVal = maxVal = matrix[0][0];
    for (i = 0; i <= userParameters->getMaxAA(); i++)
    {
        for (j = 1; j <= i; j++)
        {
            if (matrix[i][j] < minVal)
            {
                minVal = matrix[i][j];
            }
            if (matrix[i][j] > maxVal)
            {
                maxVal = matrix[i][j];
            }
        }
    }

    if (!minimise)
    {
        // Ajusta para tornar positivo, se negFlag == false
        if (!negFlag && minVal < 0)
        {
            for (i = 0; i <= userParameters->getMaxAA(); i++)
            {
                ti = (*xref)[i];
                if (ti != -1)
                {
                    for (j = 0; j <= userParameters->getMaxAA(); j++)
                    {
                        tj = (*xref)[j];
                        if (tj != -1)
                        {
                            matrix[ti][tj] -= minVal;
                        }
                    }
                }
            }
        }
        // Ajuste de gap
        int _gapPos1 = userParameters->getGapPos1();
        int _gapPos2 = userParameters->getGapPos2();

        for (i = 0; i < _gapPos1; i++)
        {
            matrix[i][_gapPos1]   = grScore;
            matrix[_gapPos1][i]   = grScore;
            matrix[i][_gapPos2]   = grScore;
            matrix[_gapPos2][i]   = grScore;
        }
        matrix[_gapPos1][_gapPos1] = ggScore;
        matrix[_gapPos2][_gapPos2] = ggScore;
        matrix[_gapPos2][_gapPos1] = ggScore;
        matrix[_gapPos1][_gapPos2] = ggScore;
    }
    else
    {
        // Minimização: matrix = maxVal - matrix
        for (i = 0; i <= userParameters->getMaxAA(); i++)
        {
            for (j = 0; j <= userParameters->getMaxAA(); j++)
            {
                matrix[i][j] = maxVal - matrix[i][j];
            }
        }
    }

    maxRes += 2;
    return maxRes;
}

bool SubMatrix::getUserMatFromFile(char *str, int alignResidueType, int alignType)
{
    checkResidueAndAlignType(alignResidueType, alignType);

    if(userParameters->getMenuFlag())
    {
        utilityObject->getStr(string("Enter name of the matrix file"), line2);
    }
    else
    {
        line2 = string(str);
    }

    if(line2.size() == 0)
    {
        return false;
    }

    FILE *infile;
    if((infile = fopen(line2.c_str(), "r")) == NULL)
    {
        utilityObject->error("Cannot find matrix file [%s]", line2.c_str());
        return false;
    }
    fclose(infile);

    strcpy(str, line2.c_str());
    mat   = getUserMatAddress(alignResidueType, alignType);
    xref  = getUserXrefAddress(alignResidueType, alignType);

    int maxRes = 0;
    if ((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        maxRes = readMatrixSeries(str, userMat, AAXref);
    }
    else
    {
        maxRes = readUserMatrix(str, *mat, *xref);
    }

    if (maxRes <= 0) return false;
    return true;
}

void SubMatrix::compareMatrices(int mat1[NUMRES][NUMRES], int mat2[NUMRES][NUMRES])
{
    int same = 1;
    for(int row = 0; row < NUMRES; row++)
    {
        for(int col = 0; col < NUMRES; col++)
        {
            if(mat1[row][col] != mat2[row][col])
            {
                same = 0;
                cout << "The row is " << row << ". The column is " << col << endl;
                break;
            }
        }
    }

    if(same == 0)
    {
        cout << "It was not the same\n";
    }
    else
    {
        cout << "It is the same\n";
    }
}

void SubMatrix::printGetMatrixResults(int mat[NUMRES][NUMRES])
{
    ofstream outfile("getmatrix.out");
    if(!outfile)
        cerr << "oops failed to open getmatrix.out !!!\n";

    for(int row = 0; row < NUMRES; row++)
    {
        for(int col = 0; col < NUMRES; col++)
        {
            if((mat[row][col] > 9) || (mat[row][col] < 0))
            {
                outfile << " " << mat[row][col] << ",";
            }
            else
            {
                outfile << "  " << mat[row][col] << ",";
            }
        }
        outfile << "\n";
    }
}

bool SubMatrix::getUserMatSeriesFromFile(char *str)
{
    int maxRes;
    FILE *infile;
    if(userParameters->getMenuFlag())
    {
        utilityObject->getStr(string("Enter name of the matrix file"), line2);
    }
    else
    {
        line2 = string(str);
    }
    if(line2.size() == 0) return false;

    if((infile = fopen(line2.c_str(), "r")) == NULL)
    {
        utilityObject->error("Cannot find matrix file [%s]", line2.c_str());
        return false;
    }
    fclose(infile);

    strcpy(str, line2.c_str());
    maxRes = readMatrixSeries(str, userMat, AAXref);
    if (maxRes <= 0) return false;
    return true;
}

int SubMatrix::readMatrixSeries(const char *fileName, Matrix& userMat, Xref& xref)
{
    if(!fileName || fileName[0] == '\0')
    {
        utilityObject->error("comparison matrix not specified");
        return 0;
    }
    FILE *fd = fopen(fileName, "r");
    if(!fd)
    {
        utilityObject->error("cannot open %s", fileName);
        return 0;
    }
    char inline1[1024];
    while (fgets(inline1, 1024, fd) != NULL)
    {
        if (commentline(inline1)) continue;
        if (utilityObject->lineType(inline1, "CLUSTAL_SERIES"))
        {
            userSeries = true;
        }
        else
        {
            userSeries = false;
        }
        break;
    }
    // se nao for series, volta e le como matriz unica
    if(!userSeries)
    {
        fclose(fd);
        return readUserMatrix(fileName, userMat, xref);
    }

    // serie
    int nmat = 0;
    matSeries.nmat = 0;
    rewind(fd); // pois ja lemos algo
    while (fgets(inline1, 1024, fd) != NULL)
    {
        if (commentline(inline1)) continue;
        if (utilityObject->lineType(inline1, "MATRIX"))
        {
            int llimit, ulimit;
            char mat_fileName[FILENAMELEN];
            if (sscanf(inline1+6, "%d %d %s", &llimit, &ulimit, mat_fileName) != 3)
            {
                utilityObject->error("Bad format in file %s\n", fileName);
                fclose(fd);
                return 0;
            }
            if (llimit < 0 || llimit > 100 || ulimit < 0 || ulimit > 100)
            {
                utilityObject->error("Bad format in file %s\n", fileName);
                fclose(fd);
                return 0;
            }
            if (ulimit <= llimit)
            {
                utilityObject->error("in file %s: lower limit > upper (%d-%d)\n",
                                     fileName, llimit, ulimit);
                fclose(fd);
                return 0;
            }
            int n = readUserMatrix(mat_fileName, userMatSeries[nmat],
                                   AAXrefseries[nmat]);
            if(n <= 0)
            {
                utilityObject->error("Bad format in matrix file %s\n", mat_fileName);
                fclose(fd);
                return 0;
            }
            matSeries.mat[nmat].llimit = llimit;
            matSeries.mat[nmat].ulimit = ulimit;
            matSeries.mat[nmat].matptr = &userMatSeries[nmat];
            matSeries.mat[nmat].AAXref = &AAXrefseries[nmat];
            nmat++;
            if(nmat >= MAXMAT)
            {
                cout << "Matrix series file has more entries than allowed (" << MAXMAT << ").\n"
                     << "Only the first " << MAXMAT << " are used.\n";
                break;
            }
        }
    }
    fclose(fd);
    matSeries.nmat = nmat;
    return (nmat > 0 ? 1 : 0);
}

int SubMatrix::readUserMatrix(const char *fileName, Matrix& userMat, Xref& xref)
{
    if (!fileName || fileName[0] == '\0')
    {
        utilityObject->error("comparison matrix not specified");
        return 0;
    }

    FILE *fd = fopen(fileName, "r");
    if(!fd)
    {
        utilityObject->error("cannot open %s", fileName);
        return 0;
    }
    char inline1[1024];
    // lê primeira linha para pegar os códigos
    int k=0;
    while (fgets(inline1, 1024, fd) != NULL)
    {
        if (commentline(inline1)) continue;
        if (utilityObject->lineType(inline1, "CLUSTAL_SERIES"))
        {
            utilityObject->error("in %s - single matrix expected.", fileName);
            fclose(fd);
            return 0;
        }
        char codes[NUMRES];
        k=0; // quantos res
        for (int j=0; j<(int)strlen(inline1); j++)
        {
            if(isalpha((int)inline1[j])) {
                codes[k++] = inline1[j];
                if (k >= NUMRES)
                {
                    utilityObject->error("too many entries in matrix %s", fileName);
                    fclose(fd);
                    return 0;
                }
            }
        }
        if(k==0)
        {
            utilityObject->error("wrong format in matrix %s", fileName);
            fclose(fd);
            return 0;
        }
        codes[k] = '\0';
        // cross reference
        for(int i=0; i<NUMRES; i++)
        {
            xref[i] = -1;
        }
        int maxRes=0;
        for(int i=0; i<k; i++)
        {
            char c1 = codes[i];
            for(int j=0; j<NUMRES; j++)
            {
                char c2 = userParameters->getAminoAcidCode(j);
                if(c1==c2)
                {
                    xref[i] = j;
                    maxRes++;
                    break;
                }
            }
            if(xref[i]==-1 && c1!='*')
            {
                utilityObject->warning("residue %c in matrix %s not recognised",
                                       c1, fileName);
            }
        }
        // agora lê pesos
        int ix=0, ix1=0;
        while (fgets(inline1, 1024, fd) != NULL)
        {
            if(inline1[0]=='\n') continue;
            if(inline1[0]=='#'||inline1[0]=='!') break;
            char *args[NUMRES+4];
            int numargs = getArgs(inline1, args, k+1);
            if(numargs<maxRes)
            {
                utilityObject->error("wrong format in matrix %s", fileName);
                fclose(fd);
                return 0;
            }
            int farg=0;
            if(isalpha(args[0][0])) farg=1;

            float scale=1.0f;
            for(unsigned a=0; a<strlen(args[farg]); a++)
            {
                if(args[farg][a]=='.')
                {
                    scale=10.0f;
                    break;
                }
            }
            for(int i=0; i<=ix; i++)
            {
                if (xref[i]!=-1)
                {
                    double f = atof(args[i+farg]);
                    userMat.push_back( static_cast<short>(f*scale) );
                    ix1++;
                }
            }
            ix++;
            if(ix>k)
            {
                // pass
            }
        }
        userMat.resize(ix1);
        fclose(fd);
        maxRes +=2;
        return maxRes; // fim
    }

    fclose(fd);
    return 0;
}

int SubMatrix::getArgs(char *inline1,char *args[],int max)
{
    char *inptr = inline1;
    int i;
    for (i=0; i<=max; i++)
    {
        args[i] = strtok(inptr, " \t\n");
        if(!args[i]) break;
        inptr=nullptr;
    }
    return i;
}

int SubMatrix::getMatrixNum()
{
    return matrixNum;
}

int SubMatrix::getDNAMatrixNum()
{
    return DNAMatrixNum;
}

int SubMatrix::getPWMatrixNum()
{
    return pwMatrixNum;
}

int SubMatrix::getPWDNAMatrixNum()
{
    return pwDNAMatrixNum;
}

void SubMatrix::setCurrentNameAndNum(string _matrixName, int _matrixNum,
                                     int alignResidueType,int alignType)
{
    checkResidueAndAlignType(alignResidueType, alignType);

    if((alignResidueType == Protein) && (alignType == Pairwise))
    {
        pwMatrixNum     = _matrixNum;
        pwMatrixName    = new string(_matrixName);
    }
    else if((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        matrixNum       = _matrixNum;
        matrixName      = new string(_matrixName);
    }
    else if((alignResidueType == DNA) && (alignType == Pairwise))
    {
        pwDNAMatrixNum  = _matrixNum;
        pwDNAMatrixName = new string(_matrixName);
    }
    else if((alignResidueType == DNA) && (alignType == MultipleAlign))
    {
        DNAMatrixNum    = _matrixNum;
        DNAMatrixName   = new string(_matrixName);
    }

#if DEBUGFULL
    if(logObject && DEBUGLOG)
    {
        ostringstream outs;
        outs << "setCurrentNameAndNum changed matrix for ("
             << alignResidueType << "," << alignType
             << ") => " << _matrixName << "\n";
        logObject->logMsg(outs.str());
    }
#endif
}

int SubMatrix::getMatrixNumForMenu(int alignResidueType, int alignType)
{
    checkResidueAndAlignType(alignResidueType, alignType);
    if((alignResidueType == Protein) && (alignType == Pairwise))
    {
        return pwMatrixNum;
    }
    else if((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        return matrixNum;
    }
    else if((alignResidueType == DNA) && (alignType == Pairwise))
    {
        return pwDNAMatrixNum;
    }
    else if((alignResidueType == DNA) && (alignType == MultipleAlign))
    {
        return DNAMatrixNum;
    }
    return -100;
}

bool SubMatrix::commentline(char* line)
{
    if (line[0] == '#') return true;
    for(int i=0; line[i] != '\n' && line[i]!=EOS; i++)
    {
        if(!isspace(line[i])) return false;
    }
    return true;
}

void SubMatrix::printInFormat(vector<short>& temp, char* name)
{
    char nameOfFile[64];
    strcpy(nameOfFile, name);
    strcat(nameOfFile, ".out");
    ofstream outfile(nameOfFile);
    if(!outfile) cerr << "ops failed to open " << nameOfFile << "\n";

    outfile << "short " << name << "[]{\n";
    int numOnCurrentLine = 1;
    int soFar=0;
    for(int i=0; i<(int)temp.size(); i++)
    {
        if(soFar==numOnCurrentLine)
        {
            outfile<<"\n";
            soFar=0;
            numOnCurrentLine++;
        }
        if((temp[i] > 9)||(temp[i]<0)) outfile<<" "<<temp[i]<<",";
        else                          outfile<<"  "<<temp[i]<<",";
        soFar++;
        if((i+1)==(int)temp.size()-1)
        {
            if((temp[i+1]>9)||(temp[i+1]<0))
                outfile<<" "<<temp[i+1]<<"};\n";
            else
                outfile<<"  "<<temp[i+1]<<"};\n";
            break;
        }
    }

    ofstream outfile2("temp.out");
    for(int i=0; i<(int)temp.size(); i++)
    {
        outfile2 << temp[i] << " ";
    }
}

void SubMatrix::printVectorToFile(vector<short>& temp, char* name)
{
    char nameOfFile[64];
    strcpy(nameOfFile, name);
    strcat(nameOfFile, ".out");
    ofstream outfile(nameOfFile);
    if(!outfile)
        cerr << "ops fail opening " << nameOfFile << "\n";

    for(int i=0; i<(int)temp.size(); i++)
    {
        if((temp[i]>9)||(temp[i]<0)) outfile<<" "<<temp[i]<<",";
        else                        outfile<<"  "<<temp[i]<<",";
    }
    outfile.close();
}

Matrix* SubMatrix::getUserMatAddress(int alignResidueType, int alignType)
{
    if((alignResidueType == Protein) && (alignType == Pairwise))
    {
        return &pwUserMat;
    }
    else if((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        return &userMat;
    }
    else if((alignResidueType == DNA) && (alignType == Pairwise))
    {
        return &pwUserDNAMat;
    }
    else if((alignResidueType == DNA) && (alignType == MultipleAlign))
    {
        return &userDNAMat;
    }
    return nullptr;
}

Xref* SubMatrix::getUserXrefAddress(int alignResidueType, int alignType)
{
    if((alignResidueType == Protein) && (alignType == Pairwise))
    {
        return &pwAAXref;
    }
    else if((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        return &AAXref;
    }
    else if((alignResidueType == DNA) && (alignType == Pairwise))
    {
        return &pwDNAXref;
    }
    else if((alignResidueType == DNA) && (alignType == MultipleAlign))
    {
        return &DNAXref;
    }
    return nullptr;
}

void SubMatrix::checkResidueAndAlignType(int alignResidueType, int alignType)
{
    if(((alignResidueType != 0) && (alignResidueType != 1)) ||
       ((alignType != 0) && (alignType != 1)))
    {
        InvalidCombination ex(alignResidueType, alignType);
        ex.whatHappened();
        exit(1);
    }
}

void SubMatrix::tempInterface(int alignResidueType, int alignType)
{
    // Exemplo simplificado
    // Esse código de teste deve ser usado com cautela
    int matrix[NUMRES][NUMRES];
    char userFile[FILENAMELEN + 1];

    userParameters->setDNAFlag(true);
    strcpy(userFile, "mat1");
    userParameters->setMenuFlag(false);

    getUserMatFromFile(userFile, DNA, Pairwise);
    setCurrentNameAndNum(userFile, 4, DNA, Pairwise);
}

int SubMatrix::getAlnScoreMatrix(int matrix[NUMRES][NUMRES])
{
    if (!this)
    {
        throw std::runtime_error("SubMatrix::getAlnScoreMatrix() - 'this' é inválido (nulo).");
    }
    // Exemplo: BLOSUM45
    int _maxNumRes = getMatrix(blosum45mtVec, &defaultAAXref, matrix, true, 100);
    return _maxNumRes;
}

void SubMatrix::getQTMatrixForHistogram(int matrix[NUMRES][NUMRES])
{
    if (!this)
    {
        throw std::runtime_error("SubMatrix::getQTMatrixForHistogram() - 'this' é inválido (nulo).");
    }
    Matrix* _matPtrLocal = nullptr;
    Xref*   _matXrefLocal = nullptr;
    int maxRes = 0;

    if(userParameters->getDNAFlag())
    {
        if (QTDNAHistMatNum == DNAUSERDEFINED)
        {
            _matPtrLocal  = &QTscoreUserDNAMatrix;
            _matXrefLocal = &QTscoreDNAXref;
        }
        else if (QTDNAHistMatNum == DNACLUSTALW)
        {
            _matPtrLocal  = clustalvdnamtVec;
            _matXrefLocal = &defaultDNAXref;
        }
        else
        {
            _matPtrLocal  = swgapdnamtVec;
            _matXrefLocal = &defaultDNAXref;
        }
    }
    else
    {
        if (QTAAHistMatNum == AAHISTIDENTITY)
        {
            _matPtrLocal  = idmatVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTAAHistMatNum == AAHISTGONNETPAM80)
        {
            _matPtrLocal  = gon80mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTAAHistMatNum == AAHISTGONNETPAM120)
        {
            _matPtrLocal  = gon120mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTAAHistMatNum == AAHISTUSER)
        {
            _matPtrLocal  = &QTscoreUserMatrix;
            _matXrefLocal = &QTscoreXref;
        }
        else if (QTAAHistMatNum == AAHISTGONNETPAM350)
        {
            _matPtrLocal  = gon350mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else
        {
            _matPtrLocal  = gon250mtVec;
            _matXrefLocal = &defaultAAXref;
        }
    }
    if(!_matPtrLocal || !_matXrefLocal)
    {
        throw std::runtime_error("getQTMatrixForHistogram() - ponteiro de matriz/xref inválido.");
    }
    maxRes = getMatrix(_matPtrLocal, _matXrefLocal, matrix, false, 100);
}

void SubMatrix::getQTMatrixForLowScoreSeg(int matrix[NUMRES][NUMRES])
{
    if (!this)
    {
        throw std::runtime_error("SubMatrix::getQTMatrixForLowScoreSeg() - 'this' é inválido (nulo).");
    }
    Matrix* _matPtrLocal = nullptr;
    Xref*   _matXrefLocal = nullptr;
    int maxRes = 0;
    int _maxAA = userParameters->getMaxAA();
    int maxiVal = 0;
    int offset;

    if(userParameters->getDNAFlag())
    {
        if (QTsegmentDNAMatNum == DNAUSERDEFINED)
        {
            _matPtrLocal  = &QTsegmentDNAMatrix;
            _matXrefLocal = &QTsegmentDNAXref;
        }
        else if (QTsegmentDNAMatNum == DNACLUSTALW)
        {
            _matPtrLocal  = clustalvdnamtVec;
            _matXrefLocal = &defaultDNAXref;
        }
        else
        {
            _matPtrLocal  = swgapdnamtVec;
            _matXrefLocal = &defaultDNAXref;
        }
        maxRes = getMatrix(_matPtrLocal, _matXrefLocal, matrix, false, 100);
        for(int i = 0; i <= _maxAA; i++)
        {
            for(int j = 0; j <= _maxAA; j++)
            {
                if(matrix[i][j] > maxiVal)
                {
                    maxiVal = matrix[i][j];
                }
            }
        }
        offset = static_cast<int>( ( static_cast<float>(maxiVal) *
                   userParameters->getQTlowScoreDNAMarkingScale() ) / 20.0f );

        for(int i = 0; i <= _maxAA; i++)
        {
            for(int j = 0; j <= _maxAA; j++)
            {
                matrix[i][j] -= offset;
            }
        }
    }
    else
    {
        if (QTsegmentAAMatNum == QTAASEGGONNETPAM80)
        {
            _matPtrLocal  = gon80mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTsegmentAAMatNum == QTAASEGGONNETPAM120)
        {
            _matPtrLocal  = gon120mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTsegmentAAMatNum == QTAASEGUSER)
        {
            _matPtrLocal  = &QTsegmentAAMatrix;
            _matXrefLocal = &QTsegmentAAXref;
        }
        else if (QTsegmentAAMatNum == QTAASEGGONNETPAM350)
        {
            _matPtrLocal  = gon350mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else
        {
            _matPtrLocal  = gon250mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        maxRes = getMatrix(_matPtrLocal, _matXrefLocal, matrix, true, 100);
    }
}

bool SubMatrix::getQTLowScoreMatFromFile(char* fileName, bool dna)
{
    if(!fileName || fileName[0]=='\0')
        return false;

    FILE *infile = fopen(fileName, "r");
    if(!infile)
    {
        utilityObject->error("Cannot find matrix file [%s]", fileName);
        return false;
    }
    fclose(infile);

    if(dna)
    {
        int maxRes = readUserMatrix(fileName, QTsegmentDNAMatrix, QTsegmentDNAXref);
        return (maxRes>0);
    }
    else
    {
        int maxRes = readUserMatrix(fileName, QTsegmentAAMatrix, QTsegmentAAXref);
        return (maxRes>0);
    }
}

bool SubMatrix::getAAScoreMatFromFile(char *str)
{
    if(!str || str[0]=='\0')
        return false;
    FILE* infile = fopen(str,"r");
    if(!infile)
    {
        utilityObject->error("Cannot find matrix file [%s]", str);
        return false;
    }
    fclose(infile);

    int maxRes = readUserMatrix(str, QTscoreUserMatrix, QTscoreXref);
    return (maxRes>0);
}

bool SubMatrix::getDNAScoreMatFromFile(char *str)
{
    if(!str || str[0]=='\0')
        return false;
    FILE* infile = fopen(str,"r");
    if(!infile)
    {
        utilityObject->error("Cannot find matrix file [%s]", str);
        return false;
    }
    fclose(infile);

    int maxRes = readUserMatrix(str, QTscoreUserDNAMatrix, QTscoreDNAXref);
    return (maxRes>0);
}

} // namespace clustalw
