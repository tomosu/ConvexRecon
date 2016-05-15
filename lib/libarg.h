/********************************************************************/
/*                                                                  */
/*                       Argument Parse Library                     */
/*                                                                  */
/********************************************************************/

#ifndef __libarg_h__
#define __libarg_h__

#include <stdio.h>
#include <stdlib.h>
#ifdef __CYGWIN32__
#include <unistd.h>
#endif
#include <string.h>

#include <string>
#include <vector>

/********************************************************************/
/*                           class def                              */
/********************************************************************/

enum { ARG_INT, ARG_UINT, ARG_REAL, ARG_STR, ARG_STDSTR, ARG_FLAG, ARG_FLAG_BOOL, ARG_HEX };

namespace libarg{
	typedef struct _Option
	{
		const char*	name;
		int			type;
		void*		val;
		int			size;
	}Option;
}

class ARGUMENT
{
	protected:
		std::vector<libarg::Option>	m_vOption;

	public:
		 ARGUMENT( void );
		~ARGUMENT( void );

		 int addOptHex( const char* option, int& val );
		 int addOptHex( const char* option, unsigned int& val );
		 int addOption( const char* option, int& val );
		 int addOption( const char* option, unsigned int& val );
		 int addOption( const char* option, double& val );
		 int addOption( const char* ootion, std::string& val );
		 int addOption( const char* option, char* val, int size );
		 int addOpFlag( const char* option, int& val );
		 int addOpFlag( const char* option, bool& val );

		 int parse( int& argc, char* argv[] );
};


int parseGeometry( const char* str, int& x, int& y, int& sx, int& sy, int& xf, int& yf );

#endif
