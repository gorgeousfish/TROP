/* Copyright (c) 2003, 2006, 2015 StataCorp LP. All rights reserved. */

/*
 * Stata Plugin Interface (SPI) initialization.
 *
 * This file implements the mandatory plugin entry point required by the
 * Stata Plugin Interface. It stores the host-provided callback table and
 * returns the plugin protocol version so Stata can verify binary compatibility.
 *
 * This is a standard SPI boilerplate file; do not modify unless the
 * Stata Plugin Interface specification changes.
 */

#include "stplugin.h"

/* Global callback table populated by Stata at plugin load time. */
ST_plugin *_stata_ ;

/*
 * pginit -- plugin entry point called by Stata.
 *
 * Stores the callback table pointer and returns the compile-time
 * plugin version constant for binary compatibility verification.
 */
STDLL pginit(ST_plugin *p)
{
	_stata_ = p ;
	return(SD_PLUGINVER) ;
}
