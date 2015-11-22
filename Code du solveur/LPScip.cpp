// LPSCIP.cpp: implementation of the CLPScip class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"



//#include "scip_exception.h"

#include <assert.h>
#include <iostream>
#include "errmes.h"
#include "LPScip.h"

#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

//////////////////////////////////////////////////////////////////////

CLPScip::CLPScip()
{
	m_rowCount = 0;
	m_rowsPreAllocated = false;

	//Creates the scip environment
	m_status = SCIPcreate(&m_env);

	//load default plugins linke, heuristics , separators, etc.
	m_status = SCIPincludeDefaultPlugins(m_env);

	m_tolerance = 0;

	m_isMip = false;

	SCIPsetMessagehdlrLogfile(m_env, "scip.log");
}

CLPScip::~CLPScip()
{
	for (int i = 0; i < m_scipVars.size(); i++) {
		SCIPreleaseVar(m_env, &m_scipVars[i]);
	}

	for (int i = 0; i < m_scipCons.size(); i++) {
		SCIPreleaseCons(m_env, &m_scipCons[i]);
	}

	m_scipVars.clear();
	m_scipCons.clear();

	m_status = SCIPfree(&m_env);
}

bool CLPScip::isLicensed()
{
	//return m_status == SCIP_OKAY;
	return true;

}

bool CLPScip::isInitialized()
{
	//return m_env != NULL;
	return true;

}

void CLPScip::setFunction(CLPFunction* function)
{
	assert(m_scipVars.size() == 0);

	m_lpFunction = function;

	// Creates an empty problem
	m_status = SCIPcreateProbBasic(m_env, "sciProblem");

	// Set the type of the problem (Min or Max)
	switch (function->getType()) {
	case lpMinFunction:
		m_status = SCIPsetObjsense(m_env, SCIP_OBJSENSE_MINIMIZE);
		break;

	case lpMaxFunction:
		m_status = SCIPsetObjsense(m_env, SCIP_OBJSENSE_MAXIMIZE);
		break;
	}

	assert(m_status == SCIP_OKAY);

	if (function->getNumCoefficients() == 0)
		return;

#if 0
	double* lbs = &function->getLowerBounds().at(0);
	double* ubs = &function->getUpperBounds().at(0);
	double* coefs = &function->getCoefficients().at(0);
	char* types = &function->getIntegers().at(0);
#endif

	for (int i = 0; i < function->getNumCoefficients(); ++i)
	{

		SCIPVar scipVar;
		SCIP_VARTYPE type;

		//Determines the type
		switch(function->getIntegers().at(i)) {
	case 'C' :
		type = SCIP_VARTYPE_CONTINUOUS;
		break;

	case 'B' :
		type = SCIP_VARTYPE_BINARY;
		break;

	case 'I' :
		type = SCIP_VARTYPE_INTEGER;
		break;

	default:
		assert(false);
		break;

#if 0
	case 'S' :
		type = SCIP_VARTYPE_IMPLINT;
		break;

	case 'N' :
		type = SCIP_VARTYPE_IMPLINT;
		break;
#endif
		}

		// Create the variable
		if (function->getVarNames().at(i) != "")
		{
			char varName[SCIP_MAXSTRLEN];
			(void) SCIPsnprintf(varName, SCIP_MAXSTRLEN - 1, function->getVarNames().at(i).c_str(), i);
			m_status = SCIPcreateVarBasic(m_env, &scipVar, varName, function->getLowerBounds().at(i), function->getUpperBounds().at(i), function->getCoefficients().at(i), type);
		}
		else
		{
			char varName[SCIP_MAXSTRLEN];
			(void) SCIPsnprintf(varName, SCIP_MAXSTRLEN - 1, "X%d", i);
			m_status = SCIPcreateVarBasic(m_env, &scipVar, varName, function->getLowerBounds().at(i), function->getUpperBounds().at(i), function->getCoefficients().at(i), type);
		}

		if((m_status = SCIPaddVar(m_env, scipVar)) != SCIP_OKAY)
		{
			assert(false);
		}
		else
		{
			m_scipVars.push_back(scipVar);
		}
	}
	/*
	else {
	if() {
	assert(false);

	}

	}
	*/
}

int CLPScip::addConstraint(CLPConstraint* c)
{

	if (!m_rowsPreAllocated) //false
	{
		if (c->getNumCoefficients() <= 0)
			return 0;

		SCIPVar *vars = new SCIPVar[c->getNumCoefficients()];
		//double coefs = new double[c->getNumCoefficients()];

		for (int i = 0; i < c->getNumCoefficients(); ++i) {
			vars[i] = m_scipVars[c->getCoefficientColumns()[i]];
		}

		double *coeffs = &c->getCoefficients().at(0);

		int index = m_rowCount;
		SCIPCons constraint = 0;

		switch (c->getType()) {
				case lpLessThanConstraint:
					m_status = SCIPcreateConsBasicLinear(m_env, &constraint, c->getName().c_str(), c->getNumCoefficients(), vars, coeffs, -(SCIPinfinity(m_env)), c->getRhs()); 
					break;

				case lpGreaterThanConstraint:
					m_status = SCIPcreateConsBasicLinear(m_env, &constraint, c->getName().c_str(), c->getNumCoefficients(), vars, coeffs, c->getRhs(), SCIPinfinity(m_env));
					break;

				case lpEqualToConstraint:
					m_status = SCIPcreateConsBasicLinear(m_env, &constraint, c->getName().c_str(), c->getNumCoefficients(), vars, coeffs, c->getRhs(), c->getRhs());
					break;
		}

		if ((m_status = SCIPaddCons(m_env, constraint)) != 1)
		{
			assert(false);
		}
		else
		{
			m_scipCons.push_back(constraint);
			m_rowCount++;
		}

		delete [] vars;
	}
	else
	{

	}
	return getStatus();
}

int CLPScip::addConstraints(
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
	const int* rowlist, //x index of coeffs
	const int* collist, //y index of coeffs
	const double* vallist //coeffs's table
	*/

	int rowindex = -1, nbVals, constraintCounter = 0;

	int *beglist = new int[rcnt];
	double *newRhs = new double[rcnt];
	char * newSense = new char[rcnt];
	char ** newRowname = new char*[rcnt];
	

	//static int grr = 0;


	//char rowName[SCIP_MAXSTRLEN];
/*
	for (int i = 0; i < numcoefs; ++i) {
		int newrowindex = rowlist[i];
		if (newrowindex == rowindex)
			continue;

		assert(newrowindex == rowindex + 1);
		rowindex = newrowindex; 
		beglist[rowindex] = i;


		TODO : continue this code to optimize
		(void) SCIPsnprintf(rowName, SCIP_MAXSTRLEN, "R%d", consCount);
		switch(sense[consCount])
		{
			case 'L' :
			SCIPcreateEmptyRow(m_env, &row, rowName, -(SCIPinfinity(m_env)), rhs[consCount], FALSE, TRUE, TRUE);
			break;

			case 'G' :
			SCIPcreateEmptyRow(m_env, &row, rowName, rhs[consCount], SCIPinfinity(m_env), FALSE, TRUE, TRUE);
			break;

			case 'E' :
			SCIPcreateEmptyRow(m_env, &row, rowName, rhs[consCount], rhs[consCount], FALSE, TRUE, TRUE);
			break;

			default:
			break;
		}

		m_scipRows.push_back(row);
		++consCount;
		


	}*/

	for (int i = 0; i < numcoefs; ++i) {
		int newrowindex = rowlist[i];
		if (newrowindex == rowindex)
			continue;

		//assert(newrowindex == rowindex + 1);

		rowindex = newrowindex;

		beglist[constraintCounter] = i;
		newRhs[constraintCounter] = rhs[rowindex];
		newSense[constraintCounter] = sense[rowindex];

		newRowname[constraintCounter] = rowname[rowindex];
		++constraintCounter;
	}

	for (int i = 0; i < constraintCounter; ++i)
	{
		if(i < constraintCounter - 1)
			nbVals = beglist[i+1] - beglist[i];
		else
			nbVals = numcoefs - beglist[i];

		if(nbVals <= 0)
			continue;


		/*int index = collist[beglist[i]];
		SCIPaddVarsToRow(m_env, m_scipRows.at(i), nbVals, &m_scipVars.at(index), &const_cast<double*>(vallist)[beglist[i]]);
		*/

		SCIPCons constraint = 0;

		SCIPVar *vars = new SCIPVar[nbVals];
		//double coefs = new double[c->getNumCoefficients()];

		for (int j = beglist[i]; j < beglist[i] + nbVals; ++j) {
			vars[j - beglist[i]] = m_scipVars.at(collist[j]);
		}

		switch(newSense[i])
		{
		case 'L' :
			m_status = SCIPcreateConsBasicLinear(m_env, &constraint, newRowname[i], nbVals, vars, &(const_cast<double*>(vallist)[beglist[i]]), -(SCIPinfinity(m_env)), newRhs[i]);
			break;
		case 'G' :
			m_status = SCIPcreateConsBasicLinear(m_env, &constraint, newRowname[i], nbVals, vars, &(const_cast<double*>(vallist)[beglist[i]]), newRhs[i], SCIPinfinity(m_env));
			break;
		case 'E' :
			m_status = SCIPcreateConsBasicLinear(m_env, &constraint, newRowname[i], nbVals, vars, &(const_cast<double*>(vallist)[beglist[i]]), newRhs[i], newRhs[i]);
			break;
		default:

			break;
		}

		if(m_status != SCIP_OKAY)
		{
			assert(false);
		}

		if ((m_status = SCIPaddCons(m_env, constraint)) != 1)
		{
			assert(false);
		}
		else
		{
			m_scipCons.push_back(constraint);
			m_rowCount++;
		}

		delete [] vars;

		
	}

	delete [] beglist;
	delete [] newRhs;
	delete [] newSense;
	delete [] newRowname;

	//solve("loooooog", "soluce.lp", "lp");
	return getStatus();
}

bool CLPScip::solve(
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


	//changes the tolerance
#if 0
	if (m_tolerance) {
		SCIPchgFeastol(m_env, m_tolerance);
		SCIPchgLpfeastol(m_env, m_tolerance, FALSE);
		SCIPchgDualfeastol(m_env, m_tolerance);
		SCIPchgBarrierconvtol(m_env, m_tolerance);

	}
#else
	if (m_tolerance) {
		SCIPsetRealParam(m_env, "numerics/feastol", m_tolerance);
		SCIPsetRealParam(m_env, "numerics/lpfeastol", m_tolerance);
		SCIPsetRealParam(m_env, "numerics/dualfeastol", m_tolerance);
		SCIPsetRealParam(m_env, "numerics/barrierconvtol", m_tolerance);
	}

#endif



	if(timeLimit >= 0)
	{
		SCIPsetRealParam(m_env, "limits/time", timeLimit);
	}

#if 1
	if(numberOfThread > 0)
	{
		SCIPsetIntParam(m_env, "lp/threads", numberOfThread);
	}
#endif

#if 1
	if(scaling == 1)
	{
		SCIPsetBoolParam(m_env, "lp/scaling", TRUE);
	}
	else
	{
		SCIPsetBoolParam(m_env, "lp/scaling", FALSE);
	}

#endif

#if 1
	switch(reSolverType) {
	case PrimalSimplexT:
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 'p');
		break;

	case DualSimplexT:
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 'd');
		break;

#if 1
	case BarrierT:
		//barrier not supported soplex yet so we use dual mode
		//m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 'b');
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 'd');
		break;
#endif

	case CrossoverPrimalT:
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 'c');
		break;

	default: //automatic simplex
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 's');
		break;
	}
	assert(getStatus() == 0);
#endif

#if 0
	
	if(doLpPresolveInInitialSolve)
	{
		
	}

#endif

#if 1
	if(problemFile != "")
	{
		writeProblem(problemFile.c_str(), NULL);
	}
#endif

	m_status = SCIPpresolve(m_env);


	m_status = SCIPsolve(m_env);
	m_solution = SCIPgetBestSol(m_env);

	//char varname[SCIP_MAXSTRLEN];
	//_snprintf(varname, SCIP_MAXSTRLEN - 1, "%s", "problemee.txt");

	//For now we don't care about args


	if (solutionFile != "") {
		writeSolution(solutionFile.c_str());
		//m_status = SCIPprintBestSol(m_env, NULL, FALSE);
		//m_status = SCIPprintStatistics(m_env, NULL);

	}

	//if(writeStatistics)
	//{
	//	SCIPprintStatistics(m_env, NULL);
	//}

	return getStatus();
}

void CLPScip::readProblem(char* pathname, char* fileType)
{
	if (pathname != "")
	{
		SCIPfree(&m_env);
		SCIPreadProb(m_env, pathname, fileType);
	}
}

void CLPScip::writeProblem(const char* pathname, const char* fileType)
{
	m_status = SCIPwriteOrigProblem(m_env, pathname, fileType, 0);

}

double CLPScip::getPrimalObjectiveValue()
{
	//return SCIPgetBestSol(m_env);
	return SCIPgetSolOrigObj(m_env, m_solution);
	//return SCIPgetPrimalbound(m_env);

}

double CLPScip::getPrimalVariableValue(int index)
{
	return SCIPgetSolVal(m_env, m_solution, m_scipVars[index]);
}


void CLPScip::writeSolution(const TCHAR *pathname)
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

#if 0
	std::cout << "Status: " << std::endl;
	writeScipError(m_status);
	writeSolutions();
#endif


}

//void CLPScip::writeSolutions()
//{
//	
//	//////////////////////////////////////////////////////////
//	// Shows Primal Solution
//	std::cout << std::endl << "primal solution:" << std::endl;
//	std::cout << "================" << std::endl << std::endl;
//
//	m_status = SCIPprintBestSol(m_env, NULL, FALSE);
//
//	// Shows Statistics
//	std::cout << std::endl << "Statistics" << std::endl;
//	std::cout << "==========" << std::endl << std::endl;
//
//	m_status = SCIPprintStatistics(m_env, NULL);
//
//}

void CLPScip::preAllocRows(int num)
{

}

void CLPScip::getRhs(int nbConstaints, double* rhs)
{
	for(int i = 0; i < nbConstaints; ++i)
	{
		*(rhs + i) = SCIPgetRhsLinear(m_env, m_scipCons.at(i));
		
	}

}

void CLPScip::changeRhs(int index, double value, bool up, bool down)
{
	if(up)
	{
		SCIPchgRhsLinear(m_env, m_scipCons.at(index), value);
	}
	else if(down)
	{
		SCIPchgLhsLinear(m_env, m_scipCons.at(index), value);
	}
	else
	{
		assert(false);
	}


}

void CLPScip::getBounds(int start, int end, double* bounds, const char* lowOrUpper)
{
#if 1
	if (stricmp(lowOrUpper, "L") == 0)
	{
		//CPXgetlb(m_env, m_lp, bounds, start, end);
		for(int i = start; i < end; ++i)
		{
			*(bounds + i - start) = SCIPgetVarLbDive(m_env, m_scipVars.at(i));
		}
	}

	else if (stricmp(lowOrUpper, "U") == 0)
	{
		for(int i = start; i < end; ++i)
		{
			*(bounds + i - start) = SCIPgetVarUbDive(m_env, m_scipVars.at(i));
		}
	}
	else
		assert(false);
#endif
}
int CLPScip::changeBds(int varIndex, const char* lowOrUpper, double bound)
{
	if (stricmp(lowOrUpper, "L") == 0)
		m_status = SCIPchgVarLb(m_env, m_scipVars.at(varIndex), bound);

	else if (stricmp(lowOrUpper, "U") == 0)
		m_status = SCIPchgVarUb(m_env, m_scipVars.at(varIndex), bound);

	else
		assert(false);
	
	return getStatus();

}

int CLPScip::changeCoefs(int rowindex, int varIndex, double coef)
{

	SCIPVar * vars = SCIPgetVarsLinear(m_env, m_scipCons.at(rowindex));
	SCIPchgVarObj(m_env, vars[varIndex], coef);
	return getStatus();

}

void CLPScip::setTolerance(double p_tolerance)
{
	m_tolerance = p_tolerance;
}

int CLPScip::getStatus()
{
	if(m_status == SCIP_OKAY)
	{	
		return 0;
	} 
	else
	{
		return 1;
	}
	//return 0;


}

bool CLPScip::solveMIP(
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

	m_isMip = 1;

	if(problemFile != "")
	{
		writeProblem(problemFile.c_str(), NULL);
	}

#if 1
	switch(reSolverType) {
	case PrimalSimplexT:
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 'p');
		break;

	case DualSimplexT:
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 'd');
		break;

	case BarrierT:
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 'b');
		break;

	case CrossoverPrimalT:
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 'c');
		break;

	default: //automatic simplex
		m_status = SCIPsetCharParam(m_env, "lp/initalgorithm", 's');
		break;
	}
	assert(getStatus() == 0);
#endif


	m_status = SCIPsolve(m_env);
	m_solution = SCIPgetBestSol(m_env);

	//char varname[SCIP_MAXSTRLEN];
	//_snprintf(varname, SCIP_MAXSTRLEN - 1, "%s", "problemee.txt");

	//For now we don't care about args


	if (solutionFile != "") {
		writeSolution(solutionFile.c_str());
		//m_status = SCIPprintBestSol(m_env, NULL, FALSE);
		//m_status = SCIPprintStatistics(m_env, NULL);

	}

	return m_status == SCIP_OKAY;


	return true;
}

bool CLPScip::solve(CLPSolveParameters params)
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

void CLPScip::setStartingSolution(std::string filename)
{


}

void CLPScip::dumpStartingSolution(std::string filename)
{


}

std::string CLPScip::getErrorString()
{
   // the following was copied from SCIPprintError
   switch(m_status)
   {
	   case SCIP_OKAY:
		  return "normal termination";
	   case SCIP_ERROR:
		  return "unspecified error";
	   case SCIP_NOMEMORY:
		  return "insufficient memory error";
	   case SCIP_READERROR:
		  return "file read error";
	   case SCIP_WRITEERROR:
		  return "file write error";
	   case SCIP_NOFILE:
		  return "file not found error";
	   case SCIP_FILECREATEERROR:
		  return "cannot create file";
	   case SCIP_LPERROR:
		  return "error in LP solver";
	   case SCIP_NOPROBLEM:
		  return "no problem exists";
	   case SCIP_INVALIDCALL:
		  return "method cannot be called at this time in solution process";
	   case SCIP_INVALIDDATA:
		  return "method cannot be called with this type of data";
	   case SCIP_INVALIDRESULT:
		  return "method returned an invalid result code";
	   case SCIP_PLUGINNOTFOUND:
		  return "a required plugin was not found";
	   case SCIP_PARAMETERUNKNOWN:
		  return "the parameter with the given name was not found";
	   case SCIP_PARAMETERWRONGTYPE:
		  return "the parameter is not of the expected type";
	   case SCIP_PARAMETERWRONGVAL:
		  return "the value is invalid for the given parameter";
	   case SCIP_KEYALREADYEXISTING:
		  return "the given key is already existing in table";
	   case SCIP_MAXDEPTHLEVEL:
		  return "maximal branching depth level exceeded";
	   default:
		  return "";
   }
}