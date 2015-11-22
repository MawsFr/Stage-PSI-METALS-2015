#ifndef LPSOLVER_H
#define LPSOLVER_H

#include "LPBase.h"
#include "lp_lib.h"

class LP_EXPORT CLPLpsolve : public CLPBase {

public :
	CLPLpsolve();
	virtual ~CLPLpsolve();

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

private:

	lprec *m_env;
	int m_status;
	double m_tolerance;

};





#endif