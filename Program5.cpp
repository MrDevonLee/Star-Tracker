//------------------------------------------------------------------------
// Name:    Claire Endres - 50%, Devon Lee - 50%
//
// Course:  CS 1430, Section 4, Fall 2020
//
// Purpose: Maintain a list of starIDs, periods, and apparent magnitude
//          and perform calculations to find and display properties
//          such as absolute magnitude, distance, hydra status, and
//          extinction distance.
//          Class Observatory is used to store data for several stars
//          in parallel arrays and perform calculations.  User input
//          includes a command, observatory number, and necessary
//          input parameters.
//          Valid commands are A (Add add a star ID along with its
//          period and apparent magnitude, R (Remove starID)
//          P (Print the list of starIDs, periods, apparent magnitude,
//          absolute magnitude, distance, hydra status, and extinction
//          distance)
//          H (Print the list of starIDs, periods, apparent magnitude,
//          absolute magnitude, distance, hydra status, and extinction
//          distance for only those stars with Hydra status)
//          S (sorts starID, apparent magnitude, and period by
//          distance from maximum to minimum.
//          StarIDs are added to list using 'A' followed by
//          an observatory #, starID, period, and apparent magnitude
//          If a starID already exists in the list, it is not added.
//          If observatory is full, a message is displayed.
//          Star ID's are removed from the list using 'R' followed
//          by the starID to be deleted. If number doeesn't exist
//          in the list, a message is displayed.
//          When command 'P' is read, calculation is performed on the
//          period and absolute magnitude to determine the apparent
//          magnitude, distance hydra status, and extinction distance.
//          All star parameters for the observatory are displayed.
//          When command 'H' is read, calculation is performed on the
//          period and absolute magnitude to determine the apparent
//          magnitude, distance hydra status, and extinction distance.
//          Hydra star parameters for the observatory are displayed.
//          When command 'S' is read, calculation is performed on the
//          period and absolute magnitude to determine the apparent
//          distance. StarID, period, and apparent magnitude are
//          sorted by distance highest to lowest.
//          Processes commands until end-of-file.
//
// Input:   List of commands, some with arguments, terminated by
//          end-of-file.
//
// Output:  A message is printed when a star ID is added to or deleted,
//          or deleted from an observatory. A table showing starID,
//          period, apparent magnitude, absolute magnitude, distance,
//          hydra status, and extinction distance is printed showing
//          either all stars or just hydra stars. A message is printed
//          when an ovservatory is sorted by star ID.
//------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
using namespace std;

const int MAX_STARS = 6;
const int STARID_WIDTH = 8;
const int PERIOD_WIDTH = 8;
const int APP_MAG_WIDTH = 11;
const int ABS_MAG_WIDTH = 11;
const int DIST_WIDTH = 12;
const int IS_HYDRA_WIDTH = 11;
const double DIST_HYDRA = 40000.00;
const double EXTINCT_COR = 0.95;
const double ABS_MAG_NO_HYDRA = 0.721;
const double ABS_MAG_HYDRA = 0.786;
const double ABS_MAG_CONST = 0.908;
const double LOG_COEF = 1.035;
const double EXP_BASE10 = 10;
const double DENOM_5 = 5;
const int HUNDRETHS_PLACE = 2;
const int ID1 = 1967;
const int ID2 = 6325;
const int ID3 = 4523;

// Allows for sequential ordering of observatories with nosequential names
enum ObsvID_Type { ID_OBSV1 = ID1, ID_OBSV2 = ID2, ID_OBSV3 = ID3 };
enum CommandType
{
   ADD_CMD = 'A', REMOVE_CMD = 'R', SORT_DISTANCE = 'S',
   PRINT_HYDRA = 'H', PRINT_CMD = 'P'
};

//Class declarations
class Observatory
{
private:
   int starID[MAX_STARS];
   double period[MAX_STARS];
   double apparentMag[MAX_STARS];
   int numStars; //Number of stars currently stored in parallel arrays

   //---------------------------------------------------------------------
   // SwapInt swaps the values in two specified index
   // positions of an arrayof integers.
   // params: (inout, in, in)
   //---------------------------------------------------------------------
   void SwapInt(int array[], int index_1, int index_2)
   {
      int value_1 = array[index_1];
      int value_2 = array[index_2];
      array[index_1] = value_2;
      array[index_2] = value_1;
      return;
   }

   //---------------------------------------------------------------------
   // SwapDouble swaps the values in two specified index positions
   // in an array of doubles
   //
   // params: (inout, in, in)
   //---------------------------------------------------------------------
   void SwapDouble(double array[], int index_1, int index_2)
   {
      double value_1 = array[index_1];
      double value_2 = array[index_2];
      array[index_1] = value_2;
      array[index_2] = value_1;
      return;
   }

   //---------------------------------------------------------------------
   // PrintHeader prints the the first three rows of the output table
   // including the title and dashed divider line
   //
   // params: ()
   //---------------------------------------------------------------------
   void PrintHeader() const
   {
      cout << setw(STARID_WIDTH) << " "
         << setw(PERIOD_WIDTH) << " "
         << setw(APP_MAG_WIDTH) << "Apparent"
         << setw(ABS_MAG_WIDTH) << "Absolute"
         << setw(DIST_WIDTH) << " "
         << setw(IS_HYDRA_WIDTH) << " "
         << setw(DIST_WIDTH) << "Extinction" << endl
         << setw(STARID_WIDTH) << "StarID"
         << setw(PERIOD_WIDTH) << "period"
         << setw(APP_MAG_WIDTH) << "Magnitude"
         << setw(ABS_MAG_WIDTH) << "Magnitude"
         << setw(DIST_WIDTH) << "Distance"
         << setw(IS_HYDRA_WIDTH) << "IsHydra"
         << setw(DIST_WIDTH) << "Distance" << endl;
      cout << "--------------------------------------------"
         << "----------------------------------" << endl;

      return;
   }

   //---------------------------------------------------------------------
   // Find finds the index of a specified star ID and returns the index
   // or -1 if the star ID is not present.
   // Params: (in)
   //---------------------------------------------------------------------
   int Find(int givenStarID) const
   {
      for (int i = 0; i < numStars; i++)
         if (starID[i] == givenStarID)
            return i;

      return -1;
   }

   //---------------------------------------------------------------------
   // AbsoluteMagnitude calculates the absolute magnitude of the
   // specified star given its period and its status as a Hydra II star.
   // Params: (in, in)
   //---------------------------------------------------------------------
   double AbsoluteMagnitude(int currentStar, bool isHydra)
   {
      // Calculates the absolute magnutide with a modified formula
      // specific to Hydra II stars if applicable
      if (isHydra == true)
         return ABS_MAG_CONST - LOG_COEF * log10(period[currentStar])
         - ABS_MAG_HYDRA;

      // Normal formula used initially for all stars; the result is kept
      // if the stars are not Hydra II stars
      return ABS_MAG_CONST - LOG_COEF * log10(period[currentStar])
         - ABS_MAG_NO_HYDRA;
   }

   //---------------------------------------------------------------------
   // Distance calculates the distance of the specified star from its
   // apparent magnitude and absolute magnitude.
   // Params: (in, in)
   //---------------------------------------------------------------------
   double Distance(int currentStar, double absMag)
   {
      return pow(EXP_BASE10, (apparentMag[currentStar] - absMag)
         / DENOM_5 + 1);
   }

   //---------------------------------------------------------------------
   // IsHydra determines if a star is in the Hydra II Galaxy by
   // determining its distance from Earth and returns a true or false
   // depending on the status of the star.
   // Params: (in)
   //---------------------------------------------------------------------
   bool IsHydra(double distance) const
   {
      if (distance > DIST_HYDRA)
         return true;

      return false;
   }

   //---------------------------------------------------------------------
   // Extinction distance calculates the extinction distance of a star
   // from its absolute distance. This is done after determining if the
   // star is in the Hydra II Galaxy.
   // Params: (in)
   //---------------------------------------------------------------------
   double ExtinctionDistance(double distance) const
   {
      return distance * EXTINCT_COR;
   }

   //---------------------------------------------------------------------
   // PrintRow prints takes the values calculated in the print function
   // and prints them with the correct widths in a row.
   //
   // params: (in, in, in, in, in)
   //---------------------------------------------------------------------
   void PrintRow(int currentStar, double absMag, double distance,
      char isHydraOut, double extinDist) const
   {
      cout << setw(STARID_WIDTH) << starID[currentStar]
         << setw(PERIOD_WIDTH) << period[currentStar]
         << setw(APP_MAG_WIDTH) << apparentMag[currentStar]
         << setw(ABS_MAG_WIDTH) << absMag
         << setw(DIST_WIDTH) << distance
         << setw(IS_HYDRA_WIDTH) << isHydraOut
         << setw(DIST_WIDTH) << extinDist << endl;
   }

   //---------------------------------------------------------------------
   // DistanceForSort performs all of the calculations necesary for the
   // sort function in one place. It is based on the code for the print
   // functions extinction distance, "Y" and "N" for Hydra stars, and
   // isHydra. It has the ability to extract more than the distance from
   // calculations.
   //
   // params: (in)
   //---------------------------------------------------------------------
   double DistanceForSort(int currentIndex)
   {
      // Assumed non-hydra to use the correct absolute magnitude
      // calculation initially
      bool isHydra = false;
      double absMag = AbsoluteMagnitude(currentIndex, isHydra);
      double distance = Distance(currentIndex, absMag);

      // If the star is a Hydra II star, the distance is recalculated
      if (IsHydra(distance) == true)
      {
         isHydra = true;
         absMag = AbsoluteMagnitude(currentIndex, isHydra);
         distance = Distance(currentIndex, absMag);
      }
      return distance;
   }

   //---------------------------------------------------------------------
   // CalculateTable performs all of the calculations necesary for the
   // PrintAll and PrintHydra functions in one place. It returns absolute
   // magnitude, distance, Hydra II Galaxy Status, and extinction distance
   // to the calling function
   //
   // params: (out,out,out,out,out,in)
   //---------------------------------------------------------------------
   void CalculateTable(double& distance_temp, double& absMag_temp,
                       char& isHydraOut, double& extinDist,
                       bool& isHydra, int currentIndex)
   {
      // Assumed non-hydra to use the correct absolute magnitude
      // calculation initially
      isHydra = false;
      isHydraOut = 'N';
      double absMag = AbsoluteMagnitude(currentIndex, isHydra);
      double distance = Distance(currentIndex, absMag);
      distance_temp = distance;
      absMag_temp = absMag;

      // If the star is a Hydra II star, the distance is recalculated
      if (IsHydra(distance_temp) == true)
      {
         isHydra = true;
         isHydraOut = 'Y';
         absMag = AbsoluteMagnitude(currentIndex, isHydra);
         distance = Distance(currentIndex, absMag);
      }
      extinDist = ExtinctionDistance(distance);
      return;
   }

public:

   //---------------------------------------------------------------------
   // Observatory initializes an observatory object with 0 stars.
   // Params: (none)
   //---------------------------------------------------------------------
   Observatory()
   {
      // Default constructor
      numStars = 0;
   }

   //---------------------------------------------------------------------
   // Add adds a star to an observatory if the star is not already
   // present and the observatory is not full. Otherwise, it indicates
   // that the star already exists or that the observatory is full.
   // Params: (none)
   //---------------------------------------------------------------------
   void Add()
   {
      int starID_in;
      double period_in;
      double apparentMag_in;
      cin >> starID_in >> period_in >> apparentMag_in;

      int index = Find(starID_in);
      if (index >= 0)
      {
         // Warning that star already exists in the parallel arrays
         cout << "StarID " << starID_in
            << " is already in the observatory." << endl;
      }
      else if (numStars >= MAX_STARS)
      {
         // Warning that the parallel arrays are full
         cout << "The list of stars is full." << endl;
      }
      else
      {
         // Adds star to parallel arrays
         starID[numStars] = starID_in;
         period[numStars] = period_in;
         apparentMag[numStars] = apparentMag_in;
         numStars++;
         cout << "StarID " << starID_in
            << " added for observatory." << endl;
      }
      return;
   }

   //---------------------------------------------------------------------
   // Remove deletes the star from the observatory if it exists and
   // shifts the remaining stars down in the list. Otherwise, it
   // indicates that the star does not exist in the observatory
   // Params: (none)
   //---------------------------------------------------------------------
   void Remove()
   {
      int starID_in;
      cin >> starID_in;

      int index = Find(starID_in);
      if (index < 0)
      {
         // Cases where star doesn not exist in the parallel arrays
         if (numStars == 0)
            cout << "There are no stars in the observatory." << endl;
         else
            cout << starID_in << " was not removed because "
            << "it is not in the observatory. " << endl;
      }
      else
      {
         // Removes star and shifts all succeeding stars back one place
         // in the arrays to keep the sort order
         for (int i = index; i < numStars; i++)
         {
            starID[i] = starID[i + 1];
            period[i] = period[i + 1];
            apparentMag[i] = apparentMag[i + 1];
         }
         numStars--;
         cout << "Removed star " << starID_in
            << " from observatory." << endl;
      }
      return;
   }

   //---------------------------------------------------------------------
   // Sort sorts arrays starID, period, and apparent magnitude by
   // distance from maximum to minimum.
   // Params: (none)
   //---------------------------------------------------------------------
   void Sort()
   {
      // Selection sort algorithm
      int maxIndex = 0;
      for (int i = 0; i < numStars - 1; i++)
      {
         maxIndex = i;
         for (int j = i + 1; j < numStars; j++)
         {
            // Calculates distances of the stars in question to determine
            // how they should be ordered relative to each other
            double currentDistance = DistanceForSort(j);
            double maxDistance = DistanceForSort(maxIndex);

            if (currentDistance > maxDistance)
               maxIndex = j;
         }

         SwapInt(starID, i, maxIndex);
         SwapDouble(period, i, maxIndex);
         SwapDouble(apparentMag, i, maxIndex);
      }

      return;
   }

   //---------------------------------------------------------------------
   // Hydra prints out all starID, period, apparent magnitude, absolute
   // magnitude, distance, if Hydra, and extinction distance only if
   // the star is a Hydra II star
   // Params: (none)
   //---------------------------------------------------------------------
   void PrintHydra()
   {
      PrintHeader();
      double distance_temp, absMag_temp, extinDist;
      char isHydraOut;
      bool isHydra;
      for (int i = 0; i < numStars; i++)
      {
         CalculateTable(distance_temp, absMag_temp,
                        isHydraOut, extinDist, isHydra, i);


         // Star is only printed if it is in the Hydra II Galaxy
         if (isHydra == true)
            PrintRow(i, absMag_temp, distance_temp,
                     isHydraOut, extinDist);
      }
      return;
   }


   //---------------------------------------------------------------------
   // PrintAll prints out all starID, period, apparent magnitude,
   // absolute magnitude, distance, if Hydra, and extinction distance
   // Params: (none)
   //---------------------------------------------------------------------
   void PrintAll()
   {
      PrintHeader();
      double distance_temp, absMag_temp, extinDist;
      char isHydraOut;
      bool isHydra;
      for (int i = 0; i < numStars; i++)
      {
         CalculateTable(distance_temp, absMag_temp,
                        isHydraOut, extinDist, isHydra, i);

         // All stars are printed
         PrintRow(i, absMag_temp, distance_temp,
                  isHydraOut, extinDist);
      }
      return;
   }
};

// Function headers
void ProcessCommandForObservatory(Observatory& obsv, CommandType command,
   ObsvID_Type ObservatoryID);

// Main
int main()
{
   // To display results in a uniform format
   cout << fixed << showpoint << setprecision(HUNDRETHS_PLACE);
   Observatory obsv1, obsv2, obsv3; // Three objects of type observatory
   char cmd;
   int id;
   cin >> cmd >> id; // Priming read
   while (cin)
   {
      // Typecasting inputted data
      CommandType command = CommandType(cmd);
      ObsvID_Type ObservatoryID = ObsvID_Type(id);

      switch (ObservatoryID) // Calls function and passes correct object
      {
         case ID_OBSV1:
            ProcessCommandForObservatory(obsv1, command, ObservatoryID);
            break;
         case ID_OBSV2:
            ProcessCommandForObservatory(obsv2, command, ObservatoryID);
            break;
         case ID_OBSV3:
            ProcessCommandForObservatory(obsv3, command, ObservatoryID);
            break;
      }
      cin >> cmd >> id; // Reads another command until end-of-file
   }
   cout << "Normal Termination of Program 5." << endl;
   return 0;
}

//------------------------------------------------------------------------
// ProcessCommandForObservatory compares the inputted command to the
// various enumerated commmand types and calls the apropriate method in
// the observatory class based on the command.
//
// params: (inout, in, in)
//------------------------------------------------------------------------
void ProcessCommandForObservatory(Observatory& obsv, CommandType command,
   ObsvID_Type ObservatoryID)
{
   int name = ObservatoryID; // Typecasts enum so that it can be printed
   switch (command)
   {
      case ADD_CMD:
         obsv.Add();
         break;
      case REMOVE_CMD:
         obsv.Remove();
         break;
      case SORT_DISTANCE:
         obsv.Sort();
         cout << "Sorted by distance for observatory "
              << name << "." << endl; // Output here to not pass IDtext
         break;
      case PRINT_HYDRA:
         cout << "Stars in the Hydra II galaxy for Observatory with ID "
              << name << ":" << endl;
         obsv.PrintHydra();
         break;
      case PRINT_CMD:
         cout << "Stars for Observatory with ID "
              << name << ":" << endl;
         obsv.PrintAll();
         break;
   }
   return;
}
