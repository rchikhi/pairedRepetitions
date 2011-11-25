#include <jni.h>
#include "sarray.h"
/*
 * Class:     SuffixArray
 * Method:    jsarray
 * Signature: (Ljava/lang/String;[II)V
 */
JNIEXPORT void JNICALL Java_SuffixArray_jsarray
	(JNIEnv *env, jclass junk, jstring s0, jintArray a0, jint n)
{
	const jbyte *s = (*env)->GetStringUTFChars(env, s0, 0);
	jint *a = (*env)->GetIntArrayElements(env, a0, 0);
	int r = bsarray(s, a, n);
	(*env)->ReleaseStringUTFChars(env, s0, s);
	(*env)->ReleaseIntArrayElements(env, a0, a, 0);
}

/*
 * Class:     SuffixArray
 * Method:    jlcp
 * Signature: (Ljava/lang/String;[I[II)V
 */
JNIEXPORT void JNICALL Java_SuffixArray_jlcp
  	(JNIEnv *env, jclass junk, jintArray a0, jstring s0, jintArray b0, jint n)
{
	const jbyte *s = (*env)->GetStringUTFChars(env, s0, 0);
	jint *a = (*env)->GetIntArrayElements(env, a0, 0);
	jint *b = (*env)->GetIntArrayElements(env, b0, 0);
	lcpa(a, s, b, n);
	(*env)->ReleaseStringUTFChars(env, s0, s);
	(*env)->ReleaseIntArrayElements(env, a0, a, 0);
	(*env)->ReleaseIntArrayElements(env, b0, b, 0);
	
}
