#
#   common.py
#
#     Copyright (C) 2012  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#     
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


# This small module exports basic, shared classes and data


################################################################################
#                               Exceptions
################################################################################

class NotImplementedException(Exception):
    def __init__(self, method, msg=None):
        self.method = method
        self.msg = msg

    def __str__(self):
        retstr = "Method {0} has not been implemented!".format(self.method)
        if msg is not None:
            retstr += "Additional message: " + self.msg
        return retstr
