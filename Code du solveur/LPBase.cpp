// LPBase.cpp: implementation of the CLPBase class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#pragma warning(disable:4786)

#include "LPBase.h"
//#include "AisParser.h"

#include <assert.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

CLPFunction::CLPFunction()
{
	setName("function");
	m_type = lpMinFunction;
}

CLPFunction::~CLPFunction()
{
	m_coefficients.clear();
	m_lowerBounds.clear();
	m_upperBounds.clear();
	m_varNames.clear();
	m_integers.clear();
}

void CLPFunction::setNumCoefficients(int size)
{
	m_coefficients.resize(size);
	m_lowerBounds.resize(size);
	m_upperBounds.resize(size);
	m_varNames.resize(size);
	m_integers.resize(size);
}

void CLPFunction::setCoefficientAt(int index, double coeff, string vname, double lb /*= 0.0F*/, double ub /*= HUGE_VAL*/, int growby /*= 10000*/, int mustBeInteger /*= 0*/)
{
#if 0 // _DEBUG
	static int s_oldIndex = -1;
	int oldIndex = s_oldIndex;
	assert(index == 0 || index == oldIndex + 1);
	s_oldIndex = index;
#endif
	if (index >= (int)m_coefficients.size()) {
		m_coefficients.resize(m_coefficients.size() + growby);
		m_lowerBounds.resize(m_coefficients.size() + growby);
		m_upperBounds.resize(m_coefficients.size() + growby);
		m_varNames.resize(m_coefficients.size() + growby);
		m_integers.resize(m_coefficients.size() + growby);
	}
	assert(index < m_coefficients.size());

	m_coefficients[index] = coeff;
	m_lowerBounds[index] = lb;
	m_upperBounds[index] = ub;
	m_varNames[index] = vname;
	if (mustBeInteger == 1)
		m_integers[index] = 'I';
	else
		m_integers[index] = 'C';
}

///////////////////////////////////////////////////////////////////////////////

CLPConstraint::CLPConstraint()
{
	m_lastget = 0;
	m_coefficientCounter = 0;
	setName("");
	m_type = lpLessThanConstraint;
	m_rhs = 0;

	m_data1 = NULL;
	m_data2 = NULL;
	m_data3 = NULL;
}

CLPConstraint::~CLPConstraint()
{
	m_coefficientColumns.clear();
	m_coefficients.clear();
}

void CLPConstraint::setNumCoefficients(int size)
{
	m_coefficientColumns.resize(size);
	m_coefficients.resize(size);
}

void CLPConstraint::setCoefficientAt(int index, double coeff, int growby /*= 10000*/)
{
	if (coeff != 0) {
	
#if 0
		if (m_coefficientCounter >= m_coefficients.size()) {
			m_coefficients.resize(m_coefficients.size() + growby);
			//m_coefficientColumns.resize(m_coefficients.size() + growby);
			m_coefficientColumns.resize(m_coefficientColumns.size() + growby);
		}
		m_coefficients[m_coefficientCounter] = coeff;
		m_coefficientColumns[m_coefficientCounter] = index;
		++m_coefficientCounter;
#else
		// Make sure m_coefficientColumns is always sorted:
		vector<int>::iterator realEnd = m_coefficientColumns.begin() + m_coefficientCounter;
		if (binary_search(m_coefficientColumns.begin(), realEnd, index)) {
			// If it already exists, we don't update the value
			return;
		}
		else {
			if (m_coefficientCounter >= (int)m_coefficients.size()) {
				m_coefficients.resize(m_coefficients.size() + growby);
				//m_coefficientColumns.resize(m_coefficients.size() + growby);
				m_coefficientColumns.resize(m_coefficientColumns.size() + growby);
			}

			int i = m_coefficientCounter;
			// Search for the correct position to keep m_coefficientColumns sorted.
			while (i > 0) {
				if (m_coefficientColumns[i - 1] <= index)
					break;

				m_coefficientColumns[i] = m_coefficientColumns[i - 1];
				m_coefficients[i] = m_coefficients[i - 1];
				--i;
			}

			m_coefficients[i] = coeff;
			m_coefficientColumns[i] = index;
			++m_coefficientCounter;
		}
#endif
	}
}

inline int compareints (const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}

double CLPConstraint::getCoefficient(int index)
{
	if (index < 0 || m_coefficientCounter == 0)
		return 0;

	if (index > m_coefficientColumns.at(m_coefficientCounter - 1))
		return 0;
#if 0
	int* res = (int*) bsearch(&index, &m_coefficientColumns[0], m_coefficientCounter, sizeof(int), compareints);
	if (res)
		return m_coefficients[*res];
	else
		return 0;
#endif

	int numColumns = (int)m_coefficientColumns.size();
	for (int i = m_lastget; i < numColumns; i++) {
		if (m_coefficientColumns.at(i) == index) {
			m_lastget = i + 1;
			return m_coefficients.at(i);
		}
		else if (m_coefficientColumns.at(i) > index) {
			m_lastget = i;
			return 0;
		}
	}
	return 0;
}

int CLPConstraint::getCoefficientColNum(int index)
{
	if (index < (int)m_coefficientColumns.size())
		return m_coefficientColumns.at(index);

	return -1;
}

bool CLPConstraint::coefficientExist(int index)
{
	unsigned int numColumns = m_coefficientColumns.size();
	for (unsigned int k = 0; k < numColumns; ++k)
		if (m_coefficientColumns[k] == index)
			return true;
	return false;
}
///////////////////////////////////////////////////////////////////////////////

void CLPSolveParametersFactory::buildInstance(ISerialObject *&instance, const MCHAR *typeName)
{
	if (strcmp(typeName, "Solver-Parameters") == 0) {
		instance = new CLPSolveParameters();
	}
	else {
		assert(false);
	}
}

CLPSolveParameters::CLPSolveParameters()
{
	m_problemType = lpT;
	m_initialSolverType = NoSimplexT;
	m_reSolverType = NoSimplexT;

	m_doLpPresolveInInitialSolve = true;
	m_doLpPresolveInReSolve = true;
	m_scaling = 1;
	m_timeLimit = -1;
	m_numberOfThread = 0;

	m_mipSolverStrategy = 1;
	m_mipIntegerTolerance = -1;
	m_mipMaxNumNodes = -1;
	m_mipMaxNumSolutions = -1;
	m_mipAllowableGap = -1;

	m_logFile = _T("");
	m_problemFile = _T("");
	m_solutionFile = _T("");

	m_nCPX_PARAM_BARCOLNZ = -1;
	m_nCPX_PARAM_BARITLIM = -1;
	m_nCPX_PARAM_BARALG = -1;
	m_nCPX_PARAM_BARSTATALG = -1;
	m_nCPX_PARAM_DEPIND = -1;
	m_nCPX_PARAM_BARORDER = -1;

	m_writeStatistics = false;
}

CLPSolveParameters::CLPSolveParameters(CLPSolveParameters& params)
{
	m_problemType = params.m_problemType;
	m_initialSolverType = params.m_initialSolverType;
	m_reSolverType = params.m_reSolverType;

	m_doLpPresolveInInitialSolve = params.m_doLpPresolveInInitialSolve;
	m_doLpPresolveInReSolve = params.m_doLpPresolveInReSolve;
	m_scaling = params.m_scaling;
	m_timeLimit = params.m_timeLimit;
	m_numberOfThread = params.m_numberOfThread;

	m_mipSolverStrategy = params.m_mipSolverStrategy;
	m_mipIntegerTolerance = params.m_mipIntegerTolerance;
	m_mipMaxNumNodes = params.m_mipMaxNumNodes;
	m_mipMaxNumSolutions = params.m_mipMaxNumSolutions;
	m_mipAllowableGap = params.m_mipAllowableGap;

	m_logFile = params.m_logFile;
	m_problemFile = params.m_problemFile;
	m_solutionFile = params.m_solutionFile;
	m_nCPX_PARAM_BARCOLNZ = params.m_nCPX_PARAM_BARCOLNZ;
	m_nCPX_PARAM_BARITLIM = params.m_nCPX_PARAM_BARITLIM;
	m_nCPX_PARAM_BARALG = params.m_nCPX_PARAM_BARALG;
	m_nCPX_PARAM_BARSTATALG = params.m_nCPX_PARAM_BARSTATALG;
	m_nCPX_PARAM_DEPIND = params.m_nCPX_PARAM_DEPIND;
	m_nCPX_PARAM_BARORDER = params.m_nCPX_PARAM_BARORDER;

	m_writeStatistics = params.m_writeStatistics;
}

CLPSolveParameters::~CLPSolveParameters()
{

}

// Serialization >>
BEGIN_ATTRIBUTES_MAP(CLPSolveParameters)
	DeclareStringAttribute(attribute_doLpPresolveInInitialSolve, "Do-Lp-Presolve-In-Initial-Solve", valuetype_bool, "", "Yes", true);
	DeclareStringAttribute(attribute_doLpPresolveInReSolve, "Do-Lp-Presolve-In-Re-Solve", valuetype_bool, "", "Yes", true);
	
	DeclareAttribute(attribute_problemFile, "Problem-File", valuetype_file, "");
	DeclareAttribute(attribute_logFile, "Log-File", valuetype_file, "");
	DeclareAttribute(attribute_solutionFile, "Solution-File", valuetype_file, "");
END_ATTRIBUTES_MAP()

ISerialObject::CSmrtPtr CLPSolveParameters::getAttributeValue(CAttributeId id)
{
	switch (id) {
	case attribute_doLpPresolveInInitialSolve:
		return wrapBool(m_doLpPresolveInInitialSolve);
		
	case attribute_doLpPresolveInReSolve:
		return wrapBool(m_doLpPresolveInReSolve);
		
	case attribute_problemFile:
		return wrapString(t2m(m_problemFile));

	case attribute_logFile:
		return wrapString(t2m(m_logFile));

	case attribute_solutionFile:
		return wrapString(t2m(m_solutionFile));

	default:
		return CSerialObjectSimple::getAttributeValue(id);
	}
}

#if 0
ISerialObjectFactory::CSmrtPtr CLPSolveParameters::getAttributeFactory(CAttributeId id, ESerialTypes type)
{
	switch (id) {
	case attribute_xxx_ptrobject:
		return ISerialObjectFactory::CSmrtPtr(new CYYYFactory());
		
	default:
		return CSerialObjectSimple::getAttributeFactory(id, type);
	}
	return CSerialObjectSimple::getAttributeFactory(id, type);
}
#endif

bool CLPSolveParameters::setAttributeValue(CAttributeId id, ISerialObject::CSmrtPtr value)
{
#if 1
	switch (id) {
	case attribute_doLpPresolveInInitialSolve:
		m_doLpPresolveInInitialSolve = convertToBool(*value);
		return true;
		
	case attribute_doLpPresolveInReSolve:
		m_doLpPresolveInReSolve = convertToBool(*value);
		return true;
		
	case attribute_problemFile:
		m_problemFile = m2t(value->str());
		return true;
			
	case attribute_logFile:
		m_logFile = m2t(value->str());
		return true;
			
	case attribute_solutionFile:
		m_solutionFile = m2t(value->str());
		return true;
			
	default:
		return CSerialObjectSimple::setAttributeValue(id, value);
	}
#else
	return CSerialObjectSimple::setAttributeValue(id, value);
#endif
}

int CLPSolveParameters::getAttributeNumPossibleValues(CAttributeId id)
{
#if 0
	switch (id) {
	case attribute_xxx_alternative:
		return static_cast<int>(num_alternatives);
		
	default:
		return 0;
	}
#else
	return 0;
#endif
}

std::string CLPSolveParameters::getAttributePossibleValue(CAttributeId id, int valueIndex)
{
#if 0
	switch (id) {
	case attribute_xxx_alternative:
		return getAlternativeName(static_cast<EAlternativeTypes>(valueIndex));
		
	default:
		return "";
	}
#else
	return "";
#endif
}

const std::mstring CLPSolveParameters::getObjectSummaryString()
{
	return "Solver Parameters";
}

// Serialization <<
#if 0
bool CLPSolveParameters::Parse(CParser* parser)
{
	if (!parser->isToken("["))
		return false;

	parser->readToken();
	if (!parser->isToken("Solver-Parameters"))
		return false;

	parser->readToken();
	while (!parser->isToken("]")) {
		if (parser->isToken("Problem-Type")) {
			parser->readToken();
			if (parser->isToken("Lp"))
				m_problemType = lpT;

			else if (parser->isToken("Mip"))
				m_problemType = mipT;

			else
				return false;
		}
		else if (parser->isToken("Initial-Solver-Type")) {
			parser->readToken();
			if (parser->isToken("No-Simplex"))
				m_initialSolverType = NoSimplexT;

			else if (parser->isToken("Primal-Simplex"))
				m_initialSolverType = PrimalSimplexT;

			else if (parser->isToken("Dual-Simplex"))
				m_initialSolverType = DualSimplexT;

			else
				return false;
		}
		else if (parser->isToken("Re-Solver-Type")) {
			parser->readToken();
			if (parser->isToken("No-Simplex"))
				m_reSolverType = NoSimplexT;

			else if (parser->isToken("Primal-Simplex"))
				m_reSolverType = PrimalSimplexT;

			else if (parser->isToken("Dual-Simplex"))
				m_reSolverType = DualSimplexT;

			else
				return false;
		}
		else if (parser->isToken("Do-Lp-Presolve-In-Initial-Solve")) {
			parser->readToken();
			if (parser->isToken("Yes"))
				m_doLpPresolveInInitialSolve = true;

			else if (parser->isToken("No"))
				m_doLpPresolveInInitialSolve = false;

			else
				return false;
		}
		else if (parser->isToken("Do-Lp-Presolve-In-Re-Solve")) {
			parser->readToken();
			if (parser->isToken("Yes"))
				m_doLpPresolveInReSolve = true;

			else if (parser->isToken("No"))
				m_doLpPresolveInReSolve = false;

			else
				return false;
		}
		else if (parser->isToken("Scaling")) {
			parser->readToken();
			m_scaling = static_cast<int>(atof(parser->getTokenBuffer()));
		}
		else if (parser->isToken("Mip-Solver-Strategy")) {
			parser->readToken();
			m_mipSolverStrategy = atoi(parser->getTokenBuffer());
		}
		else if (parser->isToken("Mip-Integer-Tolerance")) {
			parser->readToken();
			m_mipIntegerTolerance = atof(parser->getTokenBuffer());
		}
		else if (parser->isToken("Mip-Max-Num-Nodes")) {
			parser->readToken();
			m_mipMaxNumNodes = atoi(parser->getTokenBuffer());
		}
		else if (parser->isToken("Mip-Max-Num-Solutions")) {
			parser->readToken();
			m_mipMaxNumSolutions = atoi(parser->getTokenBuffer());
		}
		else if (parser->isToken("Mip-Allowable-Gap")) {
			parser->readToken();
			m_mipAllowableGap = atof(parser->getTokenBuffer());
		}
		else if (parser->isToken("Log-File")) {
			parser->readToken();
			m_logFile = m2t(parser->getTokenBuffer());
		}
		else if (parser->isToken("Problem-File")) {
			parser->readToken();
			m_problemFile = m2t(parser->getTokenBuffer());
		}
		else if (parser->isToken("Solution-File")) {
			parser->readToken();
			m_solutionFile = m2t(parser->getTokenBuffer());
		}
		else if (parser->isToken("Time-Limit")) {
			parser->readToken();
			m_timeLimit = atof(parser->getTokenBuffer());
		}

		parser->readToken();
	}

	return true;
}
#endif
	
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CLPBase::CLPBase()
{
	m_rowCount = 0;
	m_lpFunction = NULL;
	m_rowsPreAllocated = false;
}

CLPBase::~CLPBase()
{
	if (m_lpFunction)
	{
		delete m_lpFunction;
		m_lpFunction = NULL;
	}
}
