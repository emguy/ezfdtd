/* NOTE: This program is free software; you can redistribute it 
 * and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation; either 
 * version 3, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * Bugs can be reported to Yu Zhang <zhang235@mcmaster.ca>.
 *
 *     File Name : dft.h
 * Last Modified : Thu 11 Oct 2012 04:55:07 PM EDT
 */

#ifndef DFT_H
#define DFT_H


int setup_dft(char* file_name);
int update_dft(int time_index);
int write_dft(char* file_name);

#endif
