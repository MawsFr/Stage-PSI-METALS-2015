#ifndef LPSCIP_H
#define LPSCIP_H

#include "LPBase.h"

extern "C" {
struct Scip;
typedef struct Scip SCIP;

struct SCIP_Var;
typedef struct SCIP_Var SCIP_VAR;

struct SCIP_Cons;
typedef struct SCIP_Cons SCIP_CONS;

struct SCIP_Sol;
typedef struct SCIP_Sol SCIP_SOL;

struct SCIP_Row;
typedef struct SCIP_Row SCIP_ROW;
}

#include "scip/type_retcode.h"

typedef SCIP* SCIPEnv;
typedef SCIP_VAR* SCIPVar;
typedef SCIP_CONS* SCIPCons;
typedef SCIP_SOL* SCIPSol;
typedef SCIP_ROW* SCIPRow;


class LP_EXPORT CLPScip : public CLPBase {

public :
	CLPScip();
	virtual ~CLPScip();

	bool isLicensed();
	bool isInitialized();

	void setFunction(CLPFunction* function);


	int addConstraint(CLPConstraint* c);
	int addConstraints(int rcnt,
		const double* rhs,
		const char* sense,
		const double* rngval,
		char** rowname,
		int numcoefs,
		const int* rowlist,
		const int* collist,
		const double* vallist
	);

	bool solve(
		std::string logfile = "",
		std::string problemFile = "",
		std::string solutionFile = "",
		ESimplexSolverType initialSolverType = NoSimplexT,
		ESimplexSolverType reSolverType = NoSimplexT,
		bool doLpPresolveInInitialSolve = true,
		bool doLpPresolveInReSolve = true,
		int scaling = 1,
		double timeLimit = -1,
		int numberOfThread = 0,
		int nCPX_PARAM_BARCOLNZ = -1,
		int nCPX_PARAM_BARITLIM = -1,
		int nCPX_PARAM_BARALG = -1,
		int nCPX_PARAM_BARSTATALG = -1,
		int nCPX_PARAM_DEPIND = -1,
		int nCPX_PARAM_BARORDER = -1,
		bool writeStatistics = false
	);

	bool solveMIP(
		int mipSolverStrategy = 1,
		std::tstring logfile = _T(""),
		std::tstring problemFile = _T(""),
		std::tstring solutionFile = _T(""),
		ESimplexSolverType initialSolverType = NoSimplexT,
		ESimplexSolverType reSolverType = NoSimplexT,
		bool doLpPresolveInInitialSolve = true,
		bool doLpPresolveInReSolve = true,
		int scaling = 1,
		double mipIntegerTolerance = -1,
		int mipMaxNumNodes = -1,
		int mipMaxNumSolutions = -1,
		double mipAllowableGap = -1,
		double timeLimit = -1
		);

	bool solve(CLPSolveParameters params);
	void readProblem(char* pathname, char* fileType);
	void writeProblem(const char* pathname, const char* fileType);
	void writeSolution(const TCHAR* pathname);

	double getPrimalVariableValue(int index);
	double getPrimalObjectiveValue();

	void preAllocRows(int num);
	void getRhs(int nbConstaints, double* rhs);
	void changeRhs(int index, double value, bool up, bool down);
	void getBounds(int start, int end, double* bounds, const char* lowOrUpper);
	int changeBds(int varIndex, const char* lowOrUpper, double bound);
	int changeCoefs(int rowindex, int varIndex, double coef);
	void setTolerance(double tol);
	int getStatus();
	void setStartingSolution(std::string filename);
	void dumpStartingSolution(std::string filename);
	std::string getErrorString();


	// For debug
	void writeSolutions();
	void writeScipError(int m_status);
	

private :
	SCIPEnv m_env;
	SCIP_RETCODE m_status;
	std::vector<SCIPVar> m_scipVars;
	std::vector<SCIPCons> m_scipCons;
	SCIPSol m_solution;

	int m_isMip;
	double m_tolerance;
};

#endif