/*********************************************************************/
/*                                                                   */
/*                              libini                               */
/*                                                                   */
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libini.h"


/*********************************************************************/
/*                                ctor                               */
/*********************************************************************/

INICFG::INICFG( void )

	{
	entNum  = 0;
	lineNum = 0;
	_line   = 1;
	}


/*********************************************************************/
/*                                dtor                               */
/*********************************************************************/

INICFG::~INICFG( void )

	{
	}


/*********************************************************************/
/*                              AddEntry                             */
/*********************************************************************/

int INICFG::AddEntry( const char* name, int& ptr )

	{
	if ( entNum==MAX_ENTRY ) return -1;

	ent[ entNum ].name = name;
	ent[ entNum ].type = ENTRY_TYPE_INT;
	ent[ entNum ].val  = &ptr;
	ent[ entNum ].size = 4;

	entNum++;

	return 0;
	}


/*********************************************************************/
/*                              AddEntry                             */
/*********************************************************************/

int INICFG::AddEntry( const char* name, float& ptr )

	{
	if ( entNum==MAX_ENTRY ) return -1;

	ent[ entNum ].name = name;
	ent[ entNum ].type = ENTRY_TYPE_FLOAT;
	ent[ entNum ].val  = &ptr;
	ent[ entNum ].size = 4;

	entNum++;

	return 0;
	}



/*********************************************************************/
/*                              AddEntry                             */
/*********************************************************************/

int INICFG::AddEntry( const char* name, double& ptr )

	{
	if ( entNum==MAX_ENTRY ) return -1;

	ent[ entNum ].name = name;
	ent[ entNum ].type = ENTRY_TYPE_DOUBLE;
	ent[ entNum ].val  = &ptr;
	ent[ entNum ].size = 8;

	entNum++;

	return 0;
	}


/*********************************************************************/
/*                              AddEntry                             */
/*********************************************************************/

int INICFG::AddEntry( const char* name, char* ptr, int maxlen )

	{
	if ( entNum==MAX_ENTRY ) return -1;

	ent[ entNum ].name = name;
	ent[ entNum ].type = ENTRY_TYPE_STR;
	ent[ entNum ].val  = ptr;
	ent[ entNum ].size = maxlen;

	entNum++;

	return 0;
	}

/*********************************************************************/
/*                              AddEntry                             */
/*********************************************************************/

int INICFG::AddEntryIF( const char* name, int& ptr )

	{
	if ( entNum==MAX_ENTRY ) return -1;

	ent[ entNum ].name = name;
	ent[ entNum ].type = ENTRY_TYPE_INTORFLOAT;
	ent[ entNum ].val  = &ptr;
	ent[ entNum ].size = 4;

	entNum++;

	return 0;
	}

/*********************************************************************/
/*                             parse tables                          */
/*********************************************************************/

const char ISDELM[] =
	{
	/* +0 +1 +2 +3 +4 +5 +6 +7 +8 +9 +a +b +c +d +e +f  */
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	/*     !  "  #  $  %  &  '  (  )  *  +  ,  -  .  /  */
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 1, 0, 1, 1, 0,
	/*  0  1  2  3  4  5  6  7  8  9  :  ;  <  =  >  ?  */
	    1, 1, 1, 1, 1, 1, 1, 1,	1, 1, 1, 0, 0, 0, 0, 0,
	/*  @  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  */
	    0, 1, 1, 1, 1, 1, 1, 1,	1, 1, 1, 1, 1, 1, 1, 1,
	/*  P  Q  R  S  T  U  V  W  X  Y  Z  [  \  ]  ^  _  */
	    1, 1, 1, 1, 1, 1, 1, 1,	1, 1, 1, 0, 0, 0, 0, 1,
	/*  `  a  b  c  d  e  f  g  h  i  j  k  l  m  n  o  */
	    0, 1, 1, 1, 1, 1, 1, 1,	1, 1, 1, 1, 1, 1, 1, 1,
	/*  p  q  r  s  t  u  v  w  x  y  z  {  |  }  ~ 	*/
	    1, 1, 1, 1, 1, 1, 1, 1,	1, 1, 1, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	};


const char SKIPCHAR[] =
	{
	/* +0 +1 +2 +3 +4 +5 +6 +7 +8 +9 +a +b +c +d +e +f */
	    1, 0, 0, 0, 0, 0, 0, 0,	0, 1, 1, 0, 0, 1, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    1, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 1, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 1, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	    0, 0, 0, 0, 0, 0, 0, 0,	0, 0, 0, 0, 0, 0, 0, 0,
	};


/*********************************************************************/
/*                              get token                            */
/*********************************************************************/

int INICFG::getToken( FILE* fp, char* tp, int& line )

	{
	/************************************/
	/*             prologue             */
	/************************************/
	enum { MOD_NORM, MOD_COMM, MOD_CBRK, MOD_SBRK, MOD_QUOT };

	unsigned char 	buf;
	int 	retcd;
	int 	mode;
	int 	tokenLen;
	char*   pp = tp;


	/************************************/
	/*              reload              */
	/************************************/
Loop0:
	tp   = pp;
	*tp  = '\0';
	mode = MOD_NORM;
	tokenLen = 0;
	line = _line;

Loop1:
	retcd = fgetc( fp );
	if ( retcd==EOF ) return TOKEN_NG;

	buf = (unsigned char)retcd;

	if ( buf=='\n' ) _line++;


	/************************************/
	/*          find token head         */
	/************************************/
	if ( SKIPCHAR[ buf ] )
		goto Loop1;

	switch( buf )
		{
		case '#': 	mode = MOD_COMM; break;
		case '[': 	mode = MOD_SBRK; break;
		case '{': 	mode = MOD_CBRK; break;
		case '}': 	mode = MOD_CBRK; break;
		case '"': 	mode = MOD_QUOT; break;
		}

	if ( mode!=MOD_QUOT )
		{
		*tp++ = buf;
		*tp   = '\0';
		tokenLen++;
		}

	if ( mode==MOD_CBRK )
		return TOKEN_OK;


	/************************************/
	/*              reload              */
	/************************************/
	line = _line;

Loop2:
	retcd = fgetc( fp );
	if ( retcd==EOF ) return TOKEN_NG;

	buf = (unsigned char)retcd;

	if ( buf=='\n' ) _line++;


	/************************************/
	/*           copy 1 token           */
	/************************************/
	*tp++ = buf;
	*tp   = '\0';
	tokenLen++;

	if ( tokenLen>=MAX_TOKENLEN )
		return TOKEN_NG;

	switch( mode )
		{
		case MOD_COMM:
			if ( buf=='\n' )
				{
				*(tp-1) = '\0';
				goto Loop0;
				}
			break;

		case MOD_QUOT:
			if ( buf=='"' )
				{
				*(tp-1) = '\0';
				return TOKEN_OK;
				}
			break;

		case MOD_NORM:
			if ( !ISDELM[ buf ] )
				{
				*(tp-1) = '\0';
				return TOKEN_OK;
				}
			break;

		case MOD_SBRK:
			if ( buf==']' ) return TOKEN_OK;
			break;
		}

	goto Loop2;
	}


/*********************************************************************/
/*                             ParseFile                             */
/*********************************************************************/

int INICFG::ParseFile( const char* fname )

	{
	/***********************************/
	/*            prologue             */
	/***********************************/
	FILE* 	fp;
	int		retcd;
	int		lineNum = 0;

	entNumInFile = 0;

	/***********************************/
	/*            open file            */
	/***********************************/
	fp = fopen( fname, "r" );
	if ( fp==NULL ) return -1;


	/***********************************/
	/*            parse loop           */
	/***********************************/
	char token[ 256 ];
	char value[ 256 ];
	char tmp  [ 256 ];

	while( 1 )
		{
		/***********************************/
		/*          get 1st token          */
		/***********************************/
		retcd = getToken( fp, token, lineNum );

		if ( feof( fp ) )
			{
			fclose( fp );
			return 0;
			}

		if ( retcd==TOKEN_EOF )
			{
			retcd = 0;
			break;
			}

		if ( retcd==TOKEN_NG )
			{
			retcd = -1;
			break;
			}


		/***********************************/
		/*        check comment line       */
		/***********************************/
		if ( retcd=='#' )
			{
			continue;
			}


		/***********************************/
		/*            get value            */
		/***********************************/
		retcd = getToken( fp, value, lineNum );

		if ( retcd==TOKEN_EOF )
			{
			retcd = -2;
			break;
			}

		if ( retcd==TOKEN_NG )
			{
			retcd = -1;
			break;
			}


		/***********************************/
		/*            set value            */
		/***********************************/
		strcpy(tmp, token);
		int l = strlen(token) - 2;
		if (l > 0)
			{
			tmp[l] = 0;
			}
		for ( int i=0; i<entNum; i++ )
			{
			int iorf = 0;
			if ( strcmp( ent[ i ].name, token ) )
				{
				if (ent[i].type == ENTRY_TYPE_INTORFLOAT)
					{
					if (strcmp( ent[ i ].name, tmp ))
						{
						continue;
						}
					if (tmp[l+1] == 'F' || tmp[l+1] == 'f')
						{
						iorf = 1;
						}
					else if (tmp[l+1] == 'I' || tmp[l+1] == 'i')
						{
						iorf = 0;
						}
					else
						{
						continue;
						}
					}
				else
					{
					continue;
					}
				}

			entNumInFile++;

			switch( ent[ i ].type )
				{
				case ENTRY_TYPE_INT:
					*(int*)ent[ i ].val = strtoul( value, NULL, 0 );
					break;

				case ENTRY_TYPE_FLOAT:
					*(float*)ent[ i ].val = (float)strtod( value, NULL );
					break;


				case ENTRY_TYPE_DOUBLE:
					*(double*)ent[ i ].val = strtod( value, NULL );
					break;

				case ENTRY_TYPE_STR:
					strncpy( (char*)ent[ i ].val, value, ent[ i ].size );
					((char*)ent[ i ].val)[ ent[ i ].size - 1 ] = 0x0;
					break;
				case ENTRY_TYPE_INTORFLOAT:
					if (iorf == 0)
						{
						*(int*)ent[ i ].val   = strtoul( value, NULL, 0 );
						}
					else
						{
						*(float*)ent[ i ].val = (float)strtod( value, NULL );
						}
					break;
				}
			}

		}


	/***********************************/
	/*            epilogue             */
	/***********************************/
	fclose( fp );

	return -1;
	}
