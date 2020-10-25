module Results exposing (..)
import Html exposing (..)
import Html.Attributes exposing (..)


-- Each element in the search result will be a list of key value pairs
type alias SearchResults = List (List(String, String))

viewSearchResults: SearchResults -> Html msg
viewSearchResults searchResults =
   List.map (\row -> tr [] (List.map (\(key, value) -> td [][text value]) row)) searchResults
   |> table [class "table table-sm table-bordered table-responsive"]


