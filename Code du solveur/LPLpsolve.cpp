// LPSOLVER.cpp: implementation of the CLPLpsolve class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#include <assert.h>
#include <iostream>
#include "errmes.h"

#include "LPLpsolve.h"

CLPLpsolve::CLPLpsolve()
{
	m_status = 0;
	m_rowCount = 0;
}

CLPLpsolve::~CLPLpsolve()
{
	if(m_env != NULL)
	{
		delete_lp(m_env);
	}

}


bool CLPLpsolve::isLicensed()
{
	return m_status == 0;
}

bool CLPLpsolve::isInitialized()
{
	return m_env != NULL;
}

void CLPLpsolve::setFunction(CLPFunction* function)
{
	m_status = 0;

	m_lpFunction = function;
	const int nVars = function->getNumCoefficients();

	// Creates an empty problem with getNumCoefficients variables
	m_env = make_lp(0, nVars);

	if(m_env == NULL)
	{	
		m_status = 1;
	}

	//if(!set_add_rowmode(m_env, FALSE))
	//{
	//	m_status = 1;
	//}
	
	//Allowing memory for rows
	//int * colno = new int[nVars];
    REAL * row= new REAL[nVars + 1];

	//set variables names
	for (int i = 0; i < nVars; ++i)
	{
		int lpIndex = i + 1;
		row[lpIndex] = function->getCoefficients().at(i);
		//colno[i] = lpIndex;

		

		//Determines the type
		switch(function->getIntegers().at(i)) {
			case 'C' :
				m_status = set_unbounded(m_env, lpIndex);
				break;

			case 'B' :
				m_status = set_binary(m_env, lpIndex, TRUE);
				break;

			case 'I' :
				m_status = set_int(m_env, lpIndex, TRUE);
				break;

			case 'S' :
				m_status = set_semicont(m_env, lpIndex, TRUE);
				break;

			default:
				assert(false);
				break;
		}
		
		//Sets upper bound and lower bound
		set_upbo(m_env, lpIndex, function->getUpperBounds().at(i));
		set_lowbo(m_env, lpIndex, function->getLowerBounds().at(i));
		set_col_name(m_env, lpIndex, const_cast<char*>(function->getVarNames().at(i).c_str()));

	}

	//set_obj_fnex(m_env, nVars, row, colno); 
	set_obj_fn(m_env, row);

	// Set the type of the problem (Min or Max)
	switch (function->getType()) {
	case lpMinFunction:
		set_minim(m_env);
		break;

	case lpMaxFunction:
		set_maxim(m_env);
		break;
	}

	if(!set_add_rowmode(m_env, TRUE))
	{
		m_status = 1;
	}
	
	delete [] row;

}

int CLPLpsolve::addConstraint(CLPConstraint* c)
{
	if (!m_rowsPreAllocated) //false
	{
		if (c->getNumCoefficients() <= 0)
			return 0;

		REAL * row = new REAL[c->getNumCoefficients()];
		int * colno = new int[c->getNumCoefficients()];
		//double coefs = new double[c->getNumCoefficients()];

		for (int i = 0; i < c->getNumCoefficients(); ++i) {
			row[i] = c->getCoefficients().at(i);
			colno[i] = (c->getCoefficientColumns().at(i)) + 1;
		}

		//REAL *row = &c->getCoefficients().at(0);
		//int *colno = &c->getCoefficientColumns().at(0);

		char sense;
		switch (c->getType()) {
		case lpLessThanConstraint:
			sense = LE;
			break;

		case lpGreaterThanConstraint:
			sense = GE;
			break;

		case lpEqualToConstraint:
			sense = EQ;
			break;
		}
		
		
		//m_status = add_constraintex(m_env, c->getNumCoefficients(), row, colno, sense, c->getRhs());
		m_status = add_constraint(m_env, row, sense, c->getRhs());
		
		++m_rowCount;
		m_status = set_row_name(m_env, m_rowCount, const_cast<char*>(c->getName().c_str()));

		delete [] row;
		delete [] colno;
		
	}

	return getStatus();
}

int CLPLpsolve::addConstraints(
							  int rcnt,
							  const double* rhs,
							  const char* sense,
							  const double* rngval,
							  char** rowname,
							  int numcoefs,
							  const int* rowlist,
							  const int* collist,
							  const double* vallist
							  )
{

	/*
	int rcnt, //number of constraints
	const double* rhs, //constraints's bounds
	const char* sense, //constraints'sense
	const double* rngval,
	char** rowname, //constraints's names
	int numcoefs, //nonzeros variables number
	const int* rowlist, //y index of coeffs
	const int* collist, //x index of coeffs
	const double* vallist //coeffs's table
	*/

	int rowindex = -1, nbVals, constraintCounter = 0;
	int *beglist = new int[rcnt];
	int * cols = new int[numcoefs];
	double *newRhs = new double[rcnt];
	char * newSense = new char[rcnt];
	char ** newRowname = new char*[rcnt];

	for (int i = 0; i < numcoefs; ++i) {

		int newrowindex = rowlist[i];
		if (newrowindex == rowindex)
			continue;

		rowindex = newrowindex;

		beglist[constraintCounter] = i;
		newRhs[constraintCounter] = rhs[rowindex];
		newSense[constraintCounter] = sense[rowindex];

		newRowname[constraintCounter] = rowname[rowindex];
		++constraintCounter;
	}

	for(int i = 0; i < constraintCounter; ++i)
	{
		if(i < constraintCounter - 1)
			nbVals = beglist[i+1] - beglist[i];
		else
			nbVals = numcoefs - beglist[i];

		if(nbVals <= 0)
			continue;

		for (int j = beglist[i]; j < beglist[i] + nbVals; ++j) {
			cols[j - beglist[i]] = collist[j] + 1;
		}

		//m_status = set_rowex(m_env, i + 1, numcoefs, const_cast<double *>(vallist), cols);

		char senseV;
		switch (newSense[i]) {
		case 'L' :
			senseV = LE;
			break;

		case 'G' :
			senseV = GE;
			break;

		case 'E' :
			senseV = EQ;
			break;
		}
		m_status = add_constraintex(m_env, nbVals, &(const_cast<double*>(vallist)[beglist[i]]), cols, senseV ,newRhs[i]);

		//m_status = set_rh_range(m_env, i + 1, rngval[i];
		m_status = set_row_name(m_env, i + 1, newRowname[i]);
	}
	
	delete [] beglist;
	delete [] cols;
	delete [] newRhs;
	delete [] newSense;
	delete [] newRowname;


	return getStatus();
}

bool CLPLpsolve::solve(
					  std::tstring logfile /*= ""*/,
					  std::tstring problemFile /*= ""*/,
					  std::tstring solutionFile /*= ""*/,
					  ESimplexSolverType initialSolverType /*= NoSimplexT*/,
					  ESimplexSolverType reSolverType /*= NoSimplexT*/,
					  bool doLpPresolveInInitialSolve /*= true*/,
					  bool doLpPresolveInReSolve /*= true*/,
					  int scaling /*= 1*/,
					  double timeLimit /*= -1*/,
					  int numberOfThread /* = 0*/,
					  int nCPX_PARAM_BARCOLNZ /*=-1*/,
					  int nCPX_PARAM_BARITLIM /*=-1*/,
					  int nCPX_PARAM_BARALG /*=1*/,
					  int nCPX_PARAM_BARSTATALG /*=1*/,
					  int nCPX_PARAM_DEPIND /*=1*/,
					  int nCPX_PARAM_BARORDER /*=1*/,
					  bool writeStatistics /*= false */
					  )
{
	if(problemFile != "")
	{
		writeProblem(problemFile.c_str(), NULL);
	}

	set_add_rowmode(m_env, FALSE);
	//m_status = m_solver->solveProblem(m_env);
	
#if 1
	switch(reSolverType) {
	case PrimalSimplexT:
		set_simplextype(m_env, SIMPLEX_PRIMAL_PRIMAL);
		break;

	case DualSimplexT:
		set_simplextype(m_env, SIMPLEX_DUAL_DUAL);
		break;

	default: //SIMPLEX_DUAL_PRIMAL
		set_simplextype(m_env, SIMPLEX_DEFAULT);
		break;
	}
#endif

	m_status = ::solve(m_env);
	//m_status = solve_LP(m_env, NULL);

	m_status = getStatus();
	


	if (solutionFile != "") {
		writeSolution(solutionFile.c_str());
	}

	return true;
}

void CLPLpsolve::readProblem(char* pathname, char* fileType)
{
	read_LP(pathname, NORMAL, NULL);
}

void CLPLpsolve::writeProblem(const char* pathname, const char* fileType)
{
	char path[200];
	_snprintf(path, sizeof(path) - 1, pathname);
	write_lp(m_env, path);
}

double CLPLpsolve::getPrimalObjectiveValue()
{
	return get_objective(m_env);
}

double CLPLpsolve::getPrimalVariableValue(int index)
{
	return get_var_primalresult(m_env, index + 1);
}


void CLPLpsolve::writeSolution(const TCHAR *pathname)
{
	char buf[200];
	buf[sizeof(buf) - 1] = '\0';

	FILE *str = _tfopen(pathname, _T("w"));
	if (!str) {
		_snprintf(buf, sizeof(buf) - 1, "Can't create file: '%s'.", t2a(pathname).c_str());
		GlobalMessage(buf);
		return;
	}


	_snprintf(buf, sizeof(buf) - 1, "Objective value :\t%f\n\n", getPrimalObjectiveValue());
	fprintf(str, buf);
	for (int i = 0; i < m_lpFunction->getNumCoefficients(); i++) {
		_snprintf(buf, sizeof(buf) - 1, "%s : \t%f\n", m_lpFunction->getVarNames().at(i).c_str(), getPrimalVariableValue(i));
		fprintf(str, buf);
	}
	fclose(str);
}

void CLPLpsolve::preAllocRows(int num)
{

}

void CLPLpsolve::getRhs(int nbConstaints, double* rhs)
{
	*rhs = get_rh(m_env, nbConstaints + 1);
}

void CLPLpsolve::changeRhs(int index, double value, bool up, bool down)
{
		set_rh(m_env, index + 1, value);
}

void CLPLpsolve::getBounds(int start, int end, double* bounds, const char* lowOrUpper)
{
if (stricmp(lowOrUpper, "L") == 0)
	{
		//CPXgetlb(m_env, m_lp, bounds, start, end);
		for(int i = start; i < end; ++i)
		{
			*(bounds + i - start) = get_lowbo(m_env, i + 1);
		}
	}

	else if (stricmp(lowOrUpper, "U") == 0)
	{
		for(int i = start; i < end; ++i)
		{
			*(bounds + i - start) = get_upbo(m_env, i + 1);
		}
	}
	else
		assert(false);

}
int CLPLpsolve::changeBds(int varIndex, const char* lowOrUpper, double bound)
{
	if (stricmp(lowOrUpper, "L") == 0)
		m_status = set_lowbo(m_env, varIndex + 1, bound);

	else if (stricmp(lowOrUpper, "U") == 0)
		m_status = set_upbo(m_env, varIndex + 1, bound);

	else
		assert(false);

	return getStatus();
}

int CLPLpsolve::changeCoefs(int rowindex, int varIndex, double coef)
{
	set_mat(m_env, rowindex + 1, varIndex + 1, coef);
	return 1;
}

void CLPLpsolve::setTolerance(double p_tolerance)
{
	m_tolerance = p_tolerance;
}

int CLPLpsolve::getStatus()
{
	if(m_status == 1)
	{	
		return 0;
	} 
	else
	{
		return 1;
	}
}

bool CLPLpsolve::solveMIP(
						 int mipSolverStrategy/* = 1*/,
						 std::tstring logfile/* = _T("")*/,
						 std::tstring problemFile/* = _T("")*/,
						 std::tstring solutionFile/* = _T("")*/,
						 ESimplexSolverType initialSolverType/* = NoSimplexT*/,
						 ESimplexSolverType reSolverType /*= NoSimplexT*/,
						 bool doLpPresolveInInitialSolve /*= true*/,
						 bool doLpPresolveInReSolve/* = true*/,
						 int scaling/* = 1*/,
						 double mipIntegerTolerance/* = -1*/,
						 int mipMaxNumNodes/* = -1*/,
						 int mipMaxNumSolutions/* = -1*/,
						 double mipAllowableGap/* = -1*/,
						 double timeLimit /*= -1*/
						 )
{
	return true;
}

bool CLPLpsolve::solve(CLPSolveParameters params)
{

	switch (params.m_problemType) {
	case lpT:
		return solve(
			params.m_logFile,
			params.m_problemFile,
			params.m_solutionFile,
			params.m_initialSolverType,
			params.m_reSolverType,
			params.m_doLpPresolveInInitialSolve,
			params.m_doLpPresolveInReSolve,
			params.m_scaling,
			params.m_timeLimit,
			params.m_numberOfThread,
			params.m_nCPX_PARAM_BARCOLNZ,
			params.m_nCPX_PARAM_BARITLIM,
			params.m_nCPX_PARAM_BARALG,
			params.m_nCPX_PARAM_BARSTATALG,
			params.m_nCPX_PARAM_DEPIND,
			params.m_nCPX_PARAM_BARORDER,
			params.m_writeStatistics
			);

	case mipT:
		return solveMIP(
			params.m_mipSolverStrategy,
			params.m_logFile,
			params.m_problemFile,
			params.m_solutionFile,
			params.m_initialSolverType,
			params.m_reSolverType,
			params.m_doLpPresolveInInitialSolve,
			params.m_doLpPresolveInReSolve,
			params.m_scaling,
			params.m_mipIntegerTolerance,
			params.m_mipMaxNumNodes,
			params.m_mipMaxNumSolutions,
			params.m_mipAllowableGap,
			params.m_timeLimit
			);
	}
	return false;

}

void CLPLpsolve::setStartingSolution(std::string filename)
{

}

void CLPLpsolve::dumpStartingSolution(std::string filename)
{

}

std::string CLPLpsolve::getErrorString()
{
	return "";
}