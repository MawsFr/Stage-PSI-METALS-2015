// LPBase.h: interface for the CLPBase class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LPBASE_H__4DD17974_6EA0_4F8B_81DD_625AC3BDC0C7__INCLUDED_)
#define AFX_LPBASE_H__4DD17974_6EA0_4F8B_81DD_625AC3BDC0C7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "SerialObjectSimple.h"
#include "SerialObjectSimpleFactory.h"


#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#ifdef LP_EXPORTS
#define LP_EXPORT __declspec(dllexport)
#else
#define LP_EXPORT __declspec(dllimport)
#endif

class CLPSolveParameters;
class CParser;

enum EProblemType {
	lpT,
	mipT
};

enum ESimplexSolverType {
	NoSimplexT,
	PrimalSimplexT,
	DualSimplexT,
	BarrierT,
	CrossoverPrimalT,
};


class CLPSolveParametersFactory : public CSerialObjectSimpleFactory
{
public:
	void buildInstance(ISerialObject *&instance, const MCHAR *typeName = _M(""));

	std::string getObjectTypeName() { return "Solver-Parameters"; }
};


class LP_EXPORT CLPSolveParameters : public CSerialObjectSimple
{
public:
	// Constructors:
	CLPSolveParameters();
	CLPSolveParameters(CLPSolveParameters& params);
	virtual ~CLPSolveParameters();

	// Serialization >>
	std::string getObjectTypeName() { return "Solver-Parameters"; }

	enum EAttributes {
		attribute_doLpPresolveInInitialSolve,
		attribute_doLpPresolveInReSolve,

		attribute_problemFile,
		attribute_logFile,
		attribute_solutionFile,
	};

	DECLARE_ATTRIBUTES_MAP()

	ISerialObject::CSmrtPtr getAttributeValue(CAttributeId id);

	//ISerialObjectFactory::CSmrtPtr getAttributeFactory(CAttributeId id, ESerialTypes type);
	bool setAttributeValue(CAttributeId id, ISerialObject::CSmrtPtr value);

	int getAttributeNumPossibleValues(CAttributeId id);
	std::string getAttributePossibleValue(CAttributeId id, int valueIndex);

	const std::mstring getObjectSummaryString();
	// Serialization <<

	//bool Parse(CParser* parser);

	EProblemType m_problemType;
	ESimplexSolverType m_initialSolverType;
	ESimplexSolverType m_reSolverType;

	bool m_doLpPresolveInInitialSolve;
	bool m_doLpPresolveInReSolve;
	int m_scaling;
	double m_timeLimit;

	int m_mipSolverStrategy;
	double m_mipIntegerTolerance;
	int m_mipMaxNumNodes;
	int m_mipMaxNumSolutions;
	double m_mipAllowableGap;

	int m_numberOfThread;
	int m_nCPX_PARAM_BARCOLNZ;
	int m_nCPX_PARAM_BARITLIM;
	int m_nCPX_PARAM_BARALG;
	int m_nCPX_PARAM_BARSTATALG;
	int m_nCPX_PARAM_DEPIND;
	int m_nCPX_PARAM_BARORDER;

	std::tstring m_logFile;
	std::tstring m_problemFile;
	std::tstring m_solutionFile;

	bool m_writeStatistics;
};

enum lpFunctionType
{
	lpMinFunction,
	lpMaxFunction
};

// vds:
// Container for:
// - The variables of the problem:
//   - lower bound,
//   - upper bound,
//   - type (integer, real)
//   - name,
//
// TODO: >>
//   - score part
//   - slack of constraint
//   - value
// TODO: <<

// - The objective function
//   - coefficients,
//   - type (minimize, maximize)
class CLPFunction
{
public:
	LP_EXPORT CLPFunction();
	LP_EXPORT virtual ~CLPFunction();

	LP_EXPORT std::string getName() { return m_name; }
	LP_EXPORT void setName(std::string name) { m_name = name; }

	LP_EXPORT inline void setType(lpFunctionType t) { m_type = t; }
	LP_EXPORT inline lpFunctionType getType() { return m_type; }

	LP_EXPORT inline int getNumCoefficients() { return (int)m_coefficients.size(); }
	LP_EXPORT void setNumCoefficients(int size);

	LP_EXPORT void setCoefficientAt(int index, double coeff, std::string name, double lb = 0, double ub = HUGE_VAL, int growby = 500, int mustBeInteger = 0);

	LP_EXPORT inline double getCoefficient(int index) { return m_coefficients.at(index); }

	LP_EXPORT std::vector<double>& getCoefficients() { return m_coefficients; }

	LP_EXPORT std::vector<double>& getLowerBounds() { return m_lowerBounds; }
	LP_EXPORT std::vector<double>& getUpperBounds() { return m_upperBounds; }
	LP_EXPORT std::vector<std::string>& getVarNames() { return m_varNames; }
	LP_EXPORT std::vector<char>& getIntegers() { return m_integers; }

private:
	std::string m_name;
	lpFunctionType m_type;
	std::vector<double> m_coefficients;

	// Variable Name
	std::vector<std::string> m_varNames;

	// Variable Type (
	// - CONTINUOUS 'C'  continuous 
	// - BINARY     'B'  binary
	// - INTEGER    'I'  general integer
	// - SEMICONT   'S'  semi-continuous
	// - SEMIINT    'N'  semi-integer
	std::vector<char> m_integers;

	// Variable Lower and Upper Bounds:
	std::vector<double> m_lowerBounds;
	std::vector<double> m_upperBounds;
};

enum lpConstraintType
{
	lpLessThanConstraint,
	lpEqualToConstraint,
	lpGreaterThanConstraint,
};

class CLPConstraint
{
public:
	LP_EXPORT CLPConstraint();
	LP_EXPORT virtual ~CLPConstraint();

	LP_EXPORT inline std::string getName() { return m_name; }
	LP_EXPORT inline void setName(std::string name) { m_name = name; }

	LP_EXPORT inline void setType(lpConstraintType t) { m_type = t; }
	LP_EXPORT inline lpConstraintType getType() { return m_type; }

	LP_EXPORT inline int getNumCoefficients() { return m_coefficientCounter; }
	LP_EXPORT void setNumCoefficients(int size);
	LP_EXPORT void setCoefficientAt(int index, double coeff, int growby = 500);

	LP_EXPORT	double getCoefficient(int index);

	LP_EXPORT int getCoefficientColNum(int index);

	LP_EXPORT inline void setRhs(double rhs) { m_rhs = rhs; }
	LP_EXPORT inline double getRhs() { return m_rhs; }

	LP_EXPORT void setData1(void* d) { m_data1 = d;}
	LP_EXPORT void* getData1() { return m_data1; }
	
	LP_EXPORT void setData2(void* d) { m_data2 = d;}
	LP_EXPORT void* getData2() { return m_data2; }
	
	LP_EXPORT void setData3(void* d) { m_data3 = d; }
	LP_EXPORT void* getData3() { return m_data3; }
	
	LP_EXPORT void print() {
		for(int i = 0; i < (int)m_coefficientColumns.size(); ++i) 
			std::cout << m_coefficientColumns[i] << " "; 
		std::cout << "\n";
		std::cout.flush();
	}

	LP_EXPORT bool coefficientExist(int index);

	LP_EXPORT inline std::vector<double>& getCoefficients() { return m_coefficients; }
	LP_EXPORT inline std::vector<int>& getCoefficientColumns() { return m_coefficientColumns; }
	LP_EXPORT inline int& getCoefficientCounter() { return m_coefficientCounter; }

private:
	lpConstraintType m_type;
	std::vector<double> m_coefficients;
	std::vector<int> m_coefficientColumns;
	int m_coefficientCounter;
	double m_rhs;
	void* m_data1;
	void* m_data2;
	void* m_data3;
	int m_lastget;

	std::string m_name;
};

class LP_EXPORT CLPBase
{
public:
	CLPBase();
	virtual ~CLPBase();

	virtual bool isLicensed() = 0;
	virtual bool isInitialized() = 0;

	virtual void setFunction(CLPFunction* f) = NULL;

	virtual void preAllocRows(int num) = NULL;

	virtual int addConstraint(CLPConstraint* c) = NULL;
	
	virtual int addConstraints(int rcnt,
		const double* rhs,
		const char* sense,
		const double* rngval,
		char** rowname,
		int numcoefs,
		const int* rowlist,
		const int* collist,
		const double* vallist
	) = NULL;
	
	virtual void getRhs(int nbConstaints, double* rhs) = NULL;
	virtual void changeRhs(int index, double value, bool up, bool down) = NULL;

	virtual void getBounds(int start, int end, double* bounds, const char* lowOrUpper) = NULL;
	virtual	int changeBds(int varIndex, const char* lowOrUpper, double bound) = NULL;

	virtual	int changeCoefs(int rowindex, int varIndex, double coef) = NULL;

	virtual void setTolerance(double tol) = NULL;

	virtual int getStatus() = NULL; //Warning !! 0 means that there is an error and 1 no error !

	virtual double getPrimalObjectiveValue() = NULL;
	virtual double getPrimalVariableValue(int index) = NULL;

	virtual bool solve(
		std::tstring logfile = _T(""),
		std::tstring problemFile = _T(""),
		std::tstring solutionFile = _T(""),
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
		) = NULL;

	virtual bool solveMIP(
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
		) = NULL;

	virtual bool solve(CLPSolveParameters params) = NULL;

	virtual void writeSolution(const TCHAR* pathname) = NULL;
	virtual void writeProblem(const char* pathname, const char* fileType) = NULL;

	virtual void setStartingSolution(std::string filename) = NULL;
	virtual void dumpStartingSolution(std::string filename) = NULL;

	virtual std::string getErrorString() = 0;

	// Helpers:
	CLPFunction* getLPFunction() { return m_lpFunction; }

	int getRowCount() { return m_rowCount; }
	void setRowCount(int p_rowCount) { m_rowCount = p_rowCount; }

protected:
	bool m_rowsPreAllocated;
	int m_rowCount;
	CLPFunction* m_lpFunction;
};

#endif // !defined(AFX_LPBASE_H__4DD17974_6EA0_4F8B_81DD_625AC3BDC0C7__INCLUDED_)
