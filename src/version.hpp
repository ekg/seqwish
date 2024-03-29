// version.hpp: Version reflection information for seqwish builds.

// modified from https://github.com/vgteam/vg/blob/master/src/version.hpp

#include <string>
#include <unordered_map>

namespace seqwish {

	using namespace std;

/// Class for holding seqwish version.
	class Version {
	public:
		/// The Git version description of this build of seqwish
		const static string VERSION;

		/// Get only the version (like v1.7.0-68-g224e7625).
		static string get_version();

		/// Get the release Git tag version of seqwish that the current version
		/// is based on (e.g. v1.7.0-68-g224e7625 will report v1.7.0).
		static string get_release();

		/// Get the codename of our released version
		static string get_codename();

		/// Get a short one-line description of the current version of seqwish with no terminating newline.
		static string get_short();

	private:
		// Not constructable
		Version() = delete;

		/// Store all the codenames by major version
		const static unordered_map<string, string> codenames;
	};

}
