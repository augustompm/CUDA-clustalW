#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include <vector>
#include <string>
#include "../general/clustalw.h"
#include "../general/userparams.h"
#include "../general/utils.h"
#include "../general/Array2D.h"

namespace clustalw
{

using namespace std;

typedef vector<short> Xref;
typedef vector<short> Matrix;

/**
 * SubMatrix
 * Responsável por lidar com as diversas matrizes de substituição (DNA e Proteína),
 * incluindo as matrizes definidas em "matrices.h", bem como as matrizes do usuário.
 */
class SubMatrix
{
public:
    SubMatrix();
    ~SubMatrix();

    /**
     * Lê argumentos de uma linha, dividindo em tokens para leitura de matrizes.
     */
    int getArgs(char *inline1, char *args[], int max);

    /**
     * Lê (eventualmente) uma ou várias matrizes de arquivo (para DNA/Protein, Pairwise/Multiple).
     */
    bool getUserMatFromFile(char *str, int alignResidueType, int alignType);
    bool getUserMatSeriesFromFile(char *str);

    bool getAAScoreMatFromFile(char *str);
    bool getDNAScoreMatFromFile(char *str);
    bool getQTLowScoreMatFromFile(char *fileName, bool dna);

    /**
     * Ajusta a matriz (ou série) atual.
     */
    void setCurrentNameAndNum(string _matrixName, int _matrixNum,
                              int alignResidueType, int alignType);

    /**
     * Retorna o "matrixNum" atual, para o menu interativo.
     */
    int getMatrixNumForMenu(int alignResidueType, int alignType);

    /**
     * Para Pairwise Alignment (constrói a submatriz final).
     */
    int getPairwiseMatrix(int matrix[NUMRES][NUMRES], PairScaleValues &scale,
                          int &matAvg);

    /**
     * Para Multiple Alignment (constrói a submatriz dependendo de % identidade).
     */
    int getProfileAlignMatrix(int matrix[NUMRES][NUMRES], double pcid, int minLen,
                              PrfScaleValues &scaleParam, int &matAvg);

    /**
     * Matriz específica para pontuar o alinhamento final (se desejado).
     */
    int getAlnScoreMatrix(int matrix[NUMRES][NUMRES]);

    /** Acesso para o menu (MatrixNum, etc.). */
    int getMatrixNum();
    int getDNAMatrixNum();
    int getPWMatrixNum();
    int getPWDNAMatrixNum();

    /** Funções relacionadas ao ClustalQt. */
    void getQTMatrixForHistogram(int matrix[NUMRES][NUMRES]);
    int getQTAAHistMatNum() { return QTAAHistMatNum; }
    int getQTDNAHistMatNum() { return QTDNAHistMatNum; }
    void setQTAAHistMatNum(int num) { QTAAHistMatNum = num; }
    void setQTDNAHistMatNum(int num) { QTDNAHistMatNum = num; }

    void getQTMatrixForLowScoreSeg(int matrix[NUMRES][NUMRES]);
    int getQTsegmentDNAMatNum() { return QTsegmentDNAMatNum; }
    void setQTsegmentDNAMatNum(int dnaMat) { QTsegmentDNAMatNum = dnaMat; }
    int getQTsegmentAAMatNum() { return QTsegmentAAMatNum; }
    void setQTsegmentAAMatNum(int aaMat) { QTsegmentAAMatNum = aaMat; }

    /**
     * Rotina de teste.
     */
    void tempInterface(int alignResidueType, int alignType);

    /**
     * Restaura valores padrão do SubMatrix.
     */
    void setValuesToDefault();

private:
    /**
     * Monta a submatriz final usando "matPtr", salvando em "matrix[][]".
     * Se "negFlag" == false, garante que valores >=0 (offset).
     * Se "minimise" == true => "SAGA matrix" (inverte).
     */
    int getMatrix(Matrix *matPtr, Xref *xref, int matrix[NUMRES][NUMRES],
                  bool negFlag, int scale, bool minimise = false);

    /**
     * Lê várias matrizes (série) de "fileName" se houver "CLUSTAL_SERIES".
     */
    int readMatrixSeries(const char *fileName, Matrix &userMat, Xref &xref);

    /**
     * Lê uma única matriz do arquivo e carrega em "userMat", "xref".
     */
    int readUserMatrix(const char *fileName, Matrix &userMat, Xref &xref);

    /**
     * Inicializa cross references defaultAAXref e defaultDNAXref.
     */
    void setUpCrossReferences();

    /**
     * Checa se uma linha é comentário ou vazia.
     */
    bool commentline(char *line);

    /**
     * Apenas debug (imprime submatriz).
     */
    void printGetMatrixResults(int mat[NUMRES][NUMRES]);

    /**
     * Apenas debug (compara duas).
     */
    void compareMatrices(int mat1[NUMRES][NUMRES], int mat2[NUMRES][NUMRES]);

    /**
     * Imprime vector<short> no arquivo, em formato triangular.
     */
    void printInFormat(vector<short> &temp, const char *name = "tempfile.out");

    /**
     * Imprime vector<short> no arquivo (linha única).
     */
    void printVectorToFile(vector<short> &temp, const char *name = "tempfile.out");

    /**
     * Retorna ponteiro p/ userMat/pwUserMat (DNA ou Protein).
     */
    Matrix *getUserMatAddress(int alignResidueType, int alignType);

    /**
     * Retorna Xref apropriado (DNA ou Protein, pairwise ou multiple).
     */
    Xref *getUserXrefAddress(int alignResidueType, int alignType);

    /**
     * Verifica se o par (alignResidueType, alignType) é válido (0 ou 1).
     */
    void checkResidueAndAlignType(int alignResidueType, int alignType);

private:
    // se for "series" = true
    bool userSeries;

    int matrixNum;      // Prot multiple
    int DNAMatrixNum;   // DNA multiple
    int pwMatrixNum;    // Prot pairwise
    int pwDNAMatrixNum; // DNA pairwise

    // Para histogram e lowScoreSegment (ClustalQt)
    int QTAAHistMatNum;
    int QTDNAHistMatNum;
    int QTsegmentDNAMatNum;
    int QTsegmentAAMatNum;

    // nomes das matrizes
    string *matrixName;      // Protein multiple
    string *DNAMatrixName;   // DNA multiple
    string *pwMatrixName;    // Protein pairwise
    string *pwDNAMatrixName; // DNA pairwise

    // Xrefs
    Xref defaultDNAXref;
    Xref defaultAAXref;

    Xref DNAXref;     // user DNA (multiple)
    Xref AAXref;      // user Prot (multiple)
    Xref pwAAXref;    // user Prot (pairwise)
    Xref pwDNAXref;   // user DNA (pairwise)
    Xref QTscoreXref;
    Xref QTscoreDNAXref;
    Xref QTsegmentDNAXref;
    Xref QTsegmentAAXref;

    // séries de matrizes
    vector<Xref> AAXrefseries;
    UserMatrixSeries matSeries;

    // Auxiliar p/ "readUserMatrix"
    string line2;

    // “userMatSeries” se for “series”
    vector<Matrix> userMatSeries;

    // user single
    Matrix userMat;
    Matrix pwUserMat;
    Matrix userDNAMat;
    Matrix pwUserDNAMat;
    Matrix QTscoreUserMatrix;
    Matrix QTscoreUserDNAMatrix;
    Matrix QTsegmentDNAMatrix;
    Matrix QTsegmentAAMatrix;

    // Tamanho em “matrices.h”
    const int sizenAAMatrix; // 276
    const int sizeDNAMatrix; // 153

    // Vetores que armazenam as definidas em “matrices.h”
    Matrix *blosum30mtVec;
    Matrix *blosum40mtVec;
    Matrix *blosum45mtVec;
    Matrix *blosum62mt2Vec;
    Matrix *blosum80mtVec;
    Matrix *pam20mtVec;
    Matrix *pam60mtVec;
    Matrix *pam120mtVec;
    Matrix *pam350mtVec;
    Matrix *idmatVec;
    Matrix *gon40mtVec;
    Matrix *gon80mtVec;
    Matrix *gon120mtVec;
    Matrix *gon160mtVec;
    Matrix *gon250mtVec;
    Matrix *gon350mtVec;
    Matrix *clustalvdnamtVec;
    Matrix *swgapdnamtVec;

    // Cópias de ponteiro que são usadas nas funções getMatrix e etc
    Matrix *mat;
    Xref *xref;
    Matrix *_matPtr;
    Xref *_matXref;

    // média do scoring
    int matrixAvgScore;
};

} // namespace clustalw

#endif
