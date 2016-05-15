/*********************************************************************/
/*                                                                   */
/*                              libini                               */
/*                                                                   */
/*********************************************************************/

#ifndef __libini_h__
#define __libini_h__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*********************************************************************/
/*                               const                               */
/*********************************************************************/

enum { ENTRY_TYPE_INT, ENTRY_TYPE_FLOAT, ENTRY_TYPE_DOUBLE, ENTRY_TYPE_STR, ENTRY_TYPE_INTORFLOAT };

const int MAX_ENTRY 	= 1024;
const int MAX_TOKENLEN 	=  256;

const int TOKEN_EOF 	=   -1;
const int TOKEN_OK  	=    0;
const int TOKEN_NG  	=    1;


/*********************************************************************/
/*                               struct                              */
/*********************************************************************/

struct ENTRY
	{
	const char*	name;
	int			type;
	void*		val;
	int			size;
	};


/*********************************************************************/
/*                               class                               */
/*********************************************************************/

class INICFG

	{
	public:
		int entNum;

		struct ENTRY ent[ MAX_ENTRY ];

		int	lineNum;

		int	_line;

		int entNumInFile;

	public:
		 INICFG( void );
		~INICFG( void );

		int AddEntry( const char* name, int&   ptr );
		int AddEntry( const char* name, float& ptr );
	  int AddEntry( const char* name, double& ptr );
		int AddEntry( const char* name, char*  ptr, int maxlen );

		int AddEntryIF( const char* name, int&   ptr );

		int ParseFile( const char* fname );

		int getToken( FILE* fp, char* tp, int& line );

		int getEntryNum(void){return entNumInFile;}
	};


#endif
