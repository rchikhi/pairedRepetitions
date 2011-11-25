/* suffix-array package

   constructors
	SuffixArray(String s0)
	SuffixArray(String s0, int a0[])
   methods
	String toString()
	int look0(String x) return smallest index i such that 
	   s[a[j]..] < x for 0<=j<i
	int look1(String s) return smallest index i such that
	   s[a[j]..] < x or x is a prefix of s[a[j]..] for 0<=j<i
	int[] lcp() return array b such that b[i] = length of
	   the longest common prefix of s[a[i-1]..] and s[a[i]..]
   notation
	s[i..] is the suffix s.substring(i)
*/

class SuffixArray {
	public final int a[];
	public final String s;

	private static native void jsarray(String s, int a[], int n);
	private static native void jlcp(int a[], String s, int b[], int n);

	SuffixArray(String s0) {
		s = s0;
		int n = s.length();
		a = new int[n+1];
		jsarray(s, a, n);
	}

	SuffixArray(String s0, int a0[]) {
		s = s0;
		a = a0;
		for(int i=0; i<s.length(); i++)
			if(s.substring(a[i]).compareTo(
			   s.substring(a[i+1])) >= 0)
				throw new Error("bad initializer");
	}

	public String toString() {
		String r = s + " {";
		for(int i=0; i<a.length; i++) {
			if(i > 0)
				r += ",";
			r += a[i];
		}
		return r + "}";
	}

	public int[] lcp() {
		int r[] = new int[a.length];
		jlcp(a, s, r, a.length);
		return r;
	}

	public int look0(String x) {
		int bot = 0;
		int top = a.length;
		while(bot < top) {
			int i = (bot+top)/2;
			if(x.compareTo(s.substring(a[i])) > 0)
				bot = i+1;
			else
				top = i;
		}
		return bot;
	}

	public int look1(String x) {
		int bot = 0;
		int top = a.length;
		while(bot < top) {
			int i = (bot+top)/2;
			String suff = s.substring(a[i]);
			if(x.compareTo(suff)>0 || suff.startsWith(x))
				bot = i+1;
			else
				top = i;
		}
		return bot;
	}

	public static void main(String args[]) {
		final String library = 
		"/u/doug/sort/suffsort/libsarray.so";

	/* test the look routines */
		int arr[] = {5,4,2,0,3,1};
		String str = "xyxyx";
		SuffixArray sa = new SuffixArray(str, arr);
		System.out.println(sa);
		sa.testlook("a");
		sa.testlook("x");
		sa.testlook("y");
		sa.testlook("z");

	/* test the native interface */
		System.load(library);
		str = "abracadabra";
		sa = new SuffixArray(str);
		System.out.println(sa);
		sa.testlook("a");
		sa.testlook("bra");
	}

	private void testlook(String x) {
		int loc0 = look0(x);
		int loc1 = look1(x);
		System.out.println(x+" at array loc "+loc0
			+"; string locs:");
		while(loc0 < loc1)
			System.out.print(" "+a[loc0++]);
		System.out.println("");
	}
}
