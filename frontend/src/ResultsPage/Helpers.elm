module ResultsPage.Helpers exposing (..)
import ResultsPage.Types exposing (..)
import SearchPage.Types

updateSearchData : Model -> SearchPage.Types.OutMsg -> Model
updateSearchData model outMsg =
    { model
        | searchHits = Just outMsg.hits
        , searchResultRows = Just outMsg.rows
    }