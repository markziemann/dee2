port module Main exposing (..)
import Browser exposing (Document)
import Html exposing (..)
import Html.Attributes exposing (..)
import Info exposing (introduction)
import Nav exposing (navbar)
import SearchBar
import SearchBarViews exposing (..)


port consoleLog : String -> Cmd msg


type alias Model =
    { searchBar : SearchBar.Model
    }


init : ( Model, Cmd Msg )
init =
    ( { searchBar = SearchBar.init
      }
    , Cmd.none
    )



---- UPDATE ----


type Msg
    = GotSearchBarMsg SearchBar.Msg


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    case msg of
        GotSearchBarMsg message ->
            SearchBar.update message model.searchBar
            |> (\(a, b) -> ({model| searchBar = a}, Cmd.map GotSearchBarMsg b))






---- VIEW ----


view : Model -> Document Msg
view model =

    { title = "Digital Expression Explorer 2"
    , body =
         [ navbar
        , div [ class "container my-5 mx-5 mx-auto" ]
            [ Html.map GotSearchBarMsg (viewLargeSearchBar model.searchBar)

            -- Using Html.map here is suboptimal. Will be refactoring
            , Html.map GotSearchBarMsg (viewSearchResults model.searchBar.searchResults)
            ]
        , introduction
        ]
    }


subscriptions : Model -> Sub Msg
subscriptions model =
    Sub.map GotSearchBarMsg (SearchBar.subscriptions model.searchBar)



---- PROGRAM ----


main : Program () Model Msg
main =
    Browser.document
        { view = view
        , init = \_ -> init
        , update = update
        , subscriptions = subscriptions
        }
