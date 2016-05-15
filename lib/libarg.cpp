/********************************************************************/
/*                                                                  */
/*                       Argument Parse Library                     */
/*                                                                  */
/********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libarg.h"


/********************************************************************/
/*                           constructor                            */
/********************************************************************/

ARGUMENT::ARGUMENT( void )
{
}


/********************************************************************/
/*                            destructor                            */
/********************************************************************/

ARGUMENT::~ARGUMENT( void )
{
}


/********************************************************************/
/*                            addOption                             */
/********************************************************************/
int ARGUMENT::addOptHex( const char* name, int& val )
{
	libarg::Option opt;
	opt.name = name;
	opt.type = ARG_HEX;
	opt.val  = &val;
	opt.size = sizeof(val);

	m_vOption.push_back(opt);
	return 0;

}

int ARGUMENT::addOptHex( const char* name, unsigned int& val )
{
	libarg::Option opt;
	opt.name = name;
	opt.type = ARG_HEX;
	opt.val  = &val;
	opt.size = sizeof(val);

	m_vOption.push_back(opt);
	return 0;

}

int ARGUMENT::addOption( const char* name, int& val )
{
	libarg::Option opt;
	
	opt.name = name;
	opt.type = ARG_INT;
	opt.val  = &val;
	opt.size = sizeof(val);

	m_vOption.push_back(opt);

	return 0;
}

int ARGUMENT::addOption( const char* name, unsigned int& val )
{
	libarg::Option opt;
	
	opt.name = name;
	opt.type = ARG_UINT;
	opt.val  = &val;
	opt.size = sizeof(val);

	m_vOption.push_back(opt);

	return 0;
}

int ARGUMENT::addOption( const char* name, double& val )
{
	libarg::Option opt;
	
	opt.name = name;
	opt.type = ARG_REAL;
	opt.val  = &val;
	opt.size = sizeof(val);

	m_vOption.push_back(opt);

	return 0;
}

int ARGUMENT::addOption( const char* name, char* val, int size )
{
	libarg::Option opt;
	
	opt.name = name;
	opt.type = ARG_STR;
	opt.val  = val;
	opt.size = size;

	m_vOption.push_back(opt);

	return 0;
}

int ARGUMENT::addOption( const char* name, std::string& val )
{
	libarg::Option opt;
	
	opt.name = name;
	opt.type = ARG_STDSTR;
	opt.val  = &val;
	opt.size = 0;

	m_vOption.push_back(opt);
	return 0;
}

int ARGUMENT::addOpFlag( const char* name, int& val )
{

	libarg::Option opt;
	
	opt.name = name;
	opt.type = ARG_FLAG;
	opt.val  = &val;
	opt.size = sizeof(val);

	m_vOption.push_back(opt);

	return 0;
}

int ARGUMENT::addOpFlag( const char* name, bool& val )
{

	libarg::Option opt;
	
	opt.name = name;
	opt.type = ARG_FLAG_BOOL;
	opt.val  = &val;
	opt.size = sizeof(val);

	m_vOption.push_back(opt);

	return 0;
}


/********************************************************************/
/*                               parse                              */
/********************************************************************/
int ARGUMENT::parse( int& argc, char* argv[] )
{
	int i,j;
	int val;

	for ( i=1; i<argc; i++ )
	{
		std::vector<libarg::Option>::iterator it = m_vOption.begin();
		
		while(it != m_vOption.end())
		{
			if (0 != strcmp((*it).name, argv[ i ] ))
			{
				++it;
				continue;
			}
			
			switch((*it).type)
			{
				case ARG_INT:
				case ARG_UINT:
					if ( i>=argc-1 ) return -1;
					// Hex ??
					if ( *argv[ i+1 ] == 'h' || *argv[ i+1 ] == 'H' )
					{
						char* ss = argv[ i+1 ];
						val = strtoul( &ss[ 1 ], NULL, 16 );
					}
					else
					{
  						val = atoi( argv[ i+1 ] );
  					}
 					*(int*)(*it).val = val;
					argv[ i   ] = NULL;
					argv[ i+1 ] = NULL;
					i++;
					break;
				case ARG_HEX:
					if ( i>=argc-1 ) return -1;
					*(int*)(*it).val = strtoul( argv[ i+1 ], NULL, 16 );
					argv[ i   ] = NULL;
					argv[ i+1 ] = NULL;
					i++;
					break;
				case ARG_REAL:
					if ( i>=argc-1 ) return -1;
					*(double*)(*it).val = atof( argv[ i+1 ] );
					argv[ i   ] = NULL;
					argv[ i+1 ] = NULL;
					i++;
					break;
				case ARG_STR:
					{
						if ( i>=argc-1 ) return -1;
						int   size   = (*it).size;
						char* pValue = (char*)(*it).val;
					
						strncpy( pValue, argv[ i+1 ], size );
						pValue[ size-1 ] = 0x0;	// Terminater
						argv[ i   ] = NULL;
						argv[ i+1 ] = NULL;
						i++;
					}
					break;
				case ARG_STDSTR:
					{
						if ( i>=argc-1 ) return -1;
						*(std::string*)(*it).val = argv[ i+1 ];	// Terminater
						argv[ i   ] = NULL;
						argv[ i+1 ] = NULL;
						i++;
					}
					break;

				case ARG_FLAG:
					*(int*)(*it).val = 1;
					argv[ i   ] = NULL;
					break;
				case ARG_FLAG_BOOL:
					*(bool*)(*it).val = true;
					argv[ i   ] = NULL;
					break;
			}
			// Next argv
			break;
		}
	}

	// adjust argument list
	for ( j=1, i=1; i<argc; i++ )
	{
		if ( argv[ i ]==NULL ) continue;
		argv[ j++ ] = argv[ i ];
	}

	argc = j;

	return 0;
}


/********************************************************************/
/*                            parse argument                        */
/********************************************************************/
//
//  [sx]x[sy]+-[x]+-[y]
//

int parseGeometry( const char* str, int& x, int& y, int& sx, int& sy, int& xf, int& yf )

	{
	char buf[ 256 ];

	/***********************************/
	/*             parse sx            */
	/***********************************/
	buf[ 0 ] = 0x0;

	unsigned int i;

	for ( i=0; i<strlen( str ); i++ )
		{
		if ( str[ i ]=='x' ) break;

		buf[ i   ] = str[ i ];
		buf[ i+1 ] = 0x0;
		}

	if ( strlen( buf )==0 )
		{
		return -1;
		}

	sx = atoi( buf );


	/***********************************/
	/*             parse sy            */
	/***********************************/
	buf[ 0 ] = 0x0;
	i++;

	unsigned int j;

	for ( j=i; j<strlen( str ); j++ )
		{
		if ( str[ j ]=='-' ) break;
		if ( str[ j ]=='+' ) break;

		buf[ j-i   ] = str[ j ];
		buf[ j-i+1 ] = 0x0;
		}

	if ( strlen( buf )==0 )
		{
		return -1;
		}

	sy = atoi( buf );


	/***********************************/
	/*             parse xf            */
	/***********************************/
	if ( str[ j ]=='-' ) xf = -1;
	if ( str[ j ]=='+' ) xf = +1;

	j++;


	/***********************************/
	/*             parse x             */
	/***********************************/
	buf[ 0 ] = 0x0;

	for ( i=j; i<strlen( str ); i++ )
		{
		if ( str[ i ]=='-' ) break;
		if ( str[ i ]=='+' ) break;

		buf[ i-j   ] = str[ i ];
		buf[ i-j+1 ] = 0x0;
		}

	if ( strlen( buf )==0 )
		{
		return -1;
		}

	x = atoi( buf );


	/***********************************/
	/*             parse yf            */
	/***********************************/
	if ( str[ i ]=='-' ) yf = -1;
	if ( str[ i ]=='+' ) yf = +1;

	i++;


	/***********************************/
	/*             parse y             */
	/***********************************/
	buf[ 0 ] = 0x0;

	for ( j=i; j<strlen( str ); j++ )
		{
		if ( str[ j ]=='-' ) break;
		if ( str[ j ]=='+' ) break;

		buf[ j-i   ] = str[ j ];
		buf[ j-i+1 ] = 0x0;
		}

	if ( strlen( buf )==0 )
		{
		return -1;
		}

	y = atoi( buf );

	return 0;
	}


/********************************************************************/
/*                               main                               */
/********************************************************************/
#if 0
int main( int argc, char* argv[] )

	{
	int    	testint;
	double 	testreal;
	char   	teststr[ 16 ];
	int    	testflag;
	int		retcd = 0;

	ARGUMENT arg;

	arg.addOption( "-int",  testint  );
	arg.addOption( "-real", testreal );
	arg.addOption( "-str",  teststr, sizeof( teststr ) );
	arg.addOpFlag( "-flag", testflag );

	retcd = arg.parse( argc, argv );
	if ( retcd || argc!=1 )
		{
		printf( "Usage: %s\n", argv[ 0 ] );
		return -1;
		}

	printf( "testint  %d\n", testint  );
	printf( "testreal %f\n", testreal );
	printf( "teststr  %s\n", teststr  );
	printf( "testflag %d\n", testflag );

	return 0;
	}

#endif
