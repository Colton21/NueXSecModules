////////////////////////////////////////////////////////////////////////
// Class:       Helperbase
// Module Type: filter
// File:        Helperbase.h
//
////////////////////////////////////////////////////////////////////////

#ifndef HELPERBASE_H
#define HELPERBASE_H

// Basic includes that are common to almost all helpers:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

namespace lee {

class HelperBase {
public:
HelperBase() {
	_debug = false;
}
~HelperBase() = default;

/**
 * @brief      Sets the debug parameter
 *
 * @param[in]  b     Boolean, debug parameter ON (true) or OFF (false)
 */
void setDebug(bool b) {
	_debug = b;
}

private:
bool _debug;
};
}

#endif
