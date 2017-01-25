//******************************************************************************
//
// Copyright (c) 2015 Microsoft Corporation. All rights reserved.
//
// This code is licensed under the MIT License (MIT).
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//******************************************************************************

#import "Starboard/SmartTypes.h"
#import "TestFramework.h"
#import "CoreFoundation/CoreFoundation.h"
#include <CoreFoundation/CFPreferences.h>

static void testPreferencesStringToCString(const char *cBuffer, 
											CFStringEncoding encoding) {
	CFStringRef expectedString = CFStringCreateWithCString(kCFAllocatorSystemDefault, 
															cBuffer, encoding);
	char* cBufferCopy = CFPreferencesStringToCString(expectedString, encoding);
	ASSERT_TRUE(cBufferCopy != NULL);
	if (cBufferCopy != NULL) {
	    //Convert characters to string using encoding of expected path string
	    CFStringRef actualString = CFStringCreateWithCString(kCFAllocatorSystemDefault, cBufferCopy, encoding);	
		//Validate if strings have the same length
		CFIndex actualLength = CFStringGetLength(actualString);
		CFIndex expectedLength = CFStringGetLength(expectedString);
		ASSERT_EQ(expectedLength, actualLength);
		if (expectedLength == actualLength) {
			//Validate if strings are the same
			for (int i = 0; i < actualLength; i++) {
			  ASSERT_EQ(cBuffer[i], cBufferCopy[i]);
			}
		}
		ASSERT_EQ(kCFCompareEqualTo, CFStringCompare(expectedString, actualString, 0));
		free(cBufferCopy);
	}
}

//Validate if string containing plain English
//characters converts properly 
//when use different encodings
TEST(CFPreferences, ConvertEnglishUTF8StringToCString) {
    //English
	const char *cPath = "C:\\Users\\Maurice FitzGerald\\AppData\\Local";
	testPreferencesStringToCString(cPath, kCFStringEncodingUTF8);
}

TEST(CFPreferences, ConvertEnglishUTF16StringToCString) {
    //English
	const char *cPath = "C:\\Users\\Maurice FitzGerald\\AppData\\Local";
	testPreferencesStringToCString(cPath, kCFStringEncodingUTF16);
}

//Validate if string containing mix of English and Russian
//characters converts properly 
//when use different encodings
TEST(CFPreferences, ConvertEnglishRussianUTF8StringToCString) {
     //English and Russian
	const char *cPath = "C:\\Users\\Корней Чуковский\\AppData\\Local";
	testPreferencesStringToCString(cPath, kCFStringEncodingUTF8);
}

TEST(CFPreferences, ConvertEnglishRussianUTF16StringToCString) {
     //English and Russian
	const char *cPath = "C:\\Users\\Корней Чуковский\\AppData\\Local";
	testPreferencesStringToCString(cPath, kCFStringEncodingUTF16);
}

//Validate if string containing mix of English and Chinese
//characters converts properly 
//when use different encodings
TEST(CFPreferences, ConvertEnglishChineseUTF8StringToCString) {
     //English and Chinese
	const char *cPath = "C:\\Users\\厚顏無恥的企鵝\\AppData\\Local";
	testPreferencesStringToCString(cPath, kCFStringEncodingUTF8);
}

TEST(CFPreferences, ConvertEnglishChineseUTF16StringToCString) {
     //English and Chinese
	const char *cPath = "C:\\Users\\厚顏無恥的企鵝\\AppData\\Local";
	testPreferencesStringToCString(cPath, kCFStringEncodingUTF16);
}

//Validate if string containing mix of English and Arabic
//characters converts properly 
//when use different encodings
TEST(CFPreferences, ConvertEnglishArabicUTF8StringToCString) {
     //English and Arabic
	const char *cPath = "C:\\Users\\البطريق وقحة\\AppData\\Local";
	testPreferencesStringToCString(cPath, kCFStringEncodingUTF8);
}

TEST(CFPreferences, ConvertEnglishArabicUTF16StringToCString) {
     //English and Arabic
	const char *cPath = "C:\\Users\\البطريق وقحة\\AppData\\Local";
	testPreferencesStringToCString(cPath, kCFStringEncodingUTF16);
}

